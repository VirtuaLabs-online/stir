#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <obs-module.h>
#include <obs-frontend-api.h>

#include "filters/stir-midpass.h"
#include "stir-context.h"
#include "chain.h"
#include "filters/common.h"

struct biquad_stage {
    float b0, b1, b2;
    float a0, a1, a2;
    float x1, x2, y1, y2;
};

struct channel_variables {
    struct biquad_stage hp;
    struct biquad_stage lp;
    struct biquad_stage bp; /* used for biquad mode */
};

struct midpass_state {
    struct filter_base base;

    float min_cutoff, max_cutoff, q;
    float wetmix, drymix;
    int mode; /* 0 = biquad band-pass, 1 = HP+LP cascade */

    struct channel_variables *ch_state[MAX_CONTEXTS * MAX_AUDIO_CHANNELS];

    float sample_rate;
    uint32_t mask;
    size_t channels;
};

static void calculate_lowpass_stage(struct biquad_stage *s, float cutoff, float sample_rate, float q)
{
    float omega_c = 2.0f * M_PI * (cutoff / sample_rate);
    float alpha = sinf(omega_c) / (2.0f * q);
    float cos_omega_c = cosf(omega_c);

    s->b0 = (1.0f - cos_omega_c) / 2.0f;
    s->b1 = 1.0f - cos_omega_c;
    s->b2 = (1.0f - cos_omega_c) / 2.0f;

    s->a0 = 1 + alpha;
    s->a1 = -2.0f * cos_omega_c;
    s->a2 = 1.0f - alpha;

    s->b0 /= s->a0;
    s->b1 /= s->a0;
    s->b2 /= s->a0;
    s->a1 /= s->a0;
    s->a2 /= s->a0;
}

static void calculate_highpass_stage(struct biquad_stage *s, float cutoff, float sample_rate, float q)
{
    float omega_c = 2.0f * M_PI * (cutoff / sample_rate);
    float alpha = sinf(omega_c) / (2.0f * q);
    float cos_omega_c = cosf(omega_c);

    s->b0 = (1.0f + cos_omega_c) / 2.0f;
    s->b1 = -(1.0f + cos_omega_c);
    s->b2 = (1.0f + cos_omega_c) / 2.0f;

    s->a0 = 1 + alpha;
    s->a1 = -2.0f * cos_omega_c;
    s->a2 = 1.0f - alpha;

    s->b0 /= s->a0;
    s->b1 /= s->a0;
    s->b2 /= s->a0;
    s->a1 /= s->a0;
    s->a2 /= s->a0;
}

static void calculate_bandpass_stage(struct biquad_stage *s, float f0, float sample_rate, float q)
{
    float omega = 2.0f * M_PI * (f0 / sample_rate);
    float alpha = sinf(omega) / (2.0f * q);
    float cosw = cosf(omega);

    s->b0 = alpha;
    s->b1 = 0.0f;
    s->b2 = -alpha;

    s->a0 = 1.0f + alpha;
    s->a1 = -2.0f * cosw;
    s->a2 = 1.0f - alpha;

    s->b0 /= s->a0;
    s->b1 /= s->a0;
    s->b2 /= s->a0;
    s->a1 /= s->a0;
    s->a2 /= s->a0;
}

static float apply_biquad_stage(struct biquad_stage *s, float in)
{
    float wet = s->b0 * in + s->b1 * s->x1 + s->b2 * s->x2 - s->a1 * s->y1 - s->a2 * s->y2;
    s->x2 = s->x1;
    s->x1 = in;
    s->y2 = s->y1;
    s->y1 = wet;
    return wet;
}

const char *stir_midpass_get_name(void *data)
{
    UNUSED_PARAMETER(data);
    return obs_module_text("STIR Midpass");
}

void stir_midpass_destroy(void *data)
{
    struct midpass_state *state = data;
    for (size_t ch = 0; ch < MAX_CONTEXTS * MAX_AUDIO_CHANNELS; ++ch) {
        if (state->ch_state[ch])
            bfree(state->ch_state[ch]);
    }
    bfree(state);
}

void stir_midpass_update(void *data, obs_data_t *settings)
{
    struct midpass_state *state = data;
    state->wetmix = (float)obs_data_get_double(settings, "mp_wet_mix");
    state->drymix = (float)obs_data_get_double(settings, "mp_dry_mix");
    state->q = (float)obs_data_get_double(settings, "mp_q");
    state->min_cutoff = (float)obs_data_get_double(settings, "mp_min_freq");
    state->max_cutoff = (float)obs_data_get_double(settings, "mp_max_freq");
    state->mode = (int)obs_data_get_int(settings, "mp_mode");
    if (state->min_cutoff >= state->max_cutoff) {
        float t = state->min_cutoff;
        state->min_cutoff = state->max_cutoff;
        state->max_cutoff = t;
    }

    state->sample_rate = (float)audio_output_get_sample_rate(obs_get_audio());
    context_collection_t *ctx_c = stir_ctx_c_find(state->base.parent);
    if (ctx_c) {
        for (size_t c = 0; c < ctx_c->length; ++c) {
            for (size_t ch = 0; ch < state->channels; ++ch) {
                uint8_t id = stir_ctx_get_num_id(ctx_c->ctx[c]);
                const char *cid = stir_ctx_get_id(ctx_c->ctx[c]);
                size_t index = id * state->channels + ch;
                char key[32];
                snprintf(key, sizeof(key), "%s_mp_ch_%zu", cid, ch % 8u);
                if (obs_data_get_bool(settings, key)) {
                    state->mask |= (1 << index);
                    if (!state->ch_state[index]) {
                        state->ch_state[index] = bzalloc(sizeof(struct channel_variables));
                    }
                } else {
                    state->mask &= ~(1 << index);
                    if (state->ch_state[index]) {
                        bfree(state->ch_state[index]);
                        state->ch_state[index] = NULL;
                    }
                }
                if (state->ch_state[index]) {
                    if (state->mode == 0) {
                        float center = sqrtf(state->min_cutoff * state->max_cutoff);
                        float q = state->q;
                        if (q <= 0.0f) {
                            q = center / (state->max_cutoff - state->min_cutoff);
                            if (q <= 0.0f)
                                q = 0.7f;
                        }
                        calculate_bandpass_stage(&state->ch_state[index]->bp, center, state->sample_rate, q);
                    } else {
                        calculate_highpass_stage(&state->ch_state[index]->hp, state->min_cutoff, state->sample_rate, state->q);
                        calculate_lowpass_stage(&state->ch_state[index]->lp, state->max_cutoff, state->sample_rate, state->q);
                    }
                }
            }
        }
    }
}

void *stir_midpass_create(obs_data_t *settings, obs_source_t *source)
{
    UNUSED_PARAMETER(settings);
    struct midpass_state *state = bzalloc(sizeof(struct midpass_state));
    state->base.ui_id = "mp";
    state->channels = audio_output_get_channels(obs_get_audio());
    state->base.context = source;
    migrate_pre_13_config(settings, state->base.ui_id, state->base.ui_id);
    return state;
}

static void process_audio(stir_context_t *ctx, void *userdata, uint32_t samplect)
{
    struct midpass_state *state = (struct midpass_state *)userdata;
    float *buf = stir_ctx_get_buf(ctx);
    uint8_t id = stir_ctx_get_num_id(ctx);
    for (size_t i = 0; i < state->channels; ++i) {
        size_t index = id * state->channels + i;
        if (state->mask & (1 << index)) {
            struct channel_variables *channel_vars = state->ch_state[index];
            for (size_t fr = 0; fr < samplect; ++fr) {
                float in = buf[i * samplect + fr];
                float out = in;
                if (state->mode == 0) {
                    float wet = apply_biquad_stage(&channel_vars->bp, in);
                    out = (in * state->drymix) + (wet * state->wetmix);
                } else {
                    float hp_out = apply_biquad_stage(&channel_vars->hp, in);
                    float lp_out = apply_biquad_stage(&channel_vars->lp, hp_out);
                    out = (in * state->drymix) + (lp_out * state->wetmix);
                }
                buf[i * samplect + fr] = out;
            }
        }
    }
}

void stir_midpass_add(void *data, obs_source_t *source)
{
    struct midpass_state *state = data;
    state->base.parent = source;
    obs_data_t *settings = obs_source_get_settings(state->base.context);
    obs_data_t *defaults = obs_data_get_defaults(settings);
    obs_data_t *settings_safe = obs_data_create_from_json(obs_data_get_json(settings));
    obs_data_t *config = obs_data_create_from_json(obs_data_get_json(defaults));
    obs_data_apply(config, settings_safe);
    stir_midpass_update(state, config);
    obs_data_release(settings_safe);
    obs_data_release(settings);
    obs_data_release(defaults);
    obs_data_release(config);
    stir_register_filter(source, "midpass", state->base.context, process_audio, state);
}

void stir_midpass_remove(void *data, obs_source_t *source)
{
    struct midpass_state *state = data;
    stir_unregister_filter(source, state->base.context);
}

obs_properties_t *stir_midpass_properties(void *data)
{
    struct midpass_state *state = data;
    obs_properties_t *props = obs_properties_create();

    filter_make_ctx_dropdown(props, &state->base);
    filter_make_ch_list(props, &state->base);

    obs_property_t *minf = obs_properties_add_float_slider(props, "mp_min_freq", "Min Frequency", 10.0, 5000.0, 1.0);
    obs_property_float_set_suffix(minf, " Hz");
    obs_property_t *maxf = obs_properties_add_float_slider(props, "mp_max_freq", "Max Frequency", 100.0, 20000.0, 1.0);
    obs_property_float_set_suffix(maxf, " Hz");
    obs_property_t *mode = obs_properties_add_list(props, "mp_mode", "Mode", OBS_COMBO_TYPE_LIST, OBS_COMBO_FORMAT_INT);
    obs_property_list_add_int(mode, "Biquad band-pass", 0);
    obs_property_list_add_int(mode, "HP + LP cascade", 1);
    obs_property_t *q = obs_properties_add_float_slider(props, "mp_q", "Q", 0.1, 10.0, 0.01);
    obs_property_float_set_suffix(q, "x");
    obs_property_t *wm = obs_properties_add_float_slider(props, "mp_wet_mix", "Wet Mix", 0.0, 1.0, 0.01);
    obs_property_float_set_suffix(wm, "x");
    obs_property_t *dm = obs_properties_add_float_slider(props, "mp_dry_mix", "Dry Mix", 0.0, 1.0, 0.01);
    obs_property_float_set_suffix(dm, "x");
    return props;
}

void stir_midpass_defaults(obs_data_t *settings)
{
    obs_data_set_default_double(settings, "mp_min_freq", 100.0);
    obs_data_set_default_double(settings, "mp_max_freq", 2000.0);
    obs_data_set_default_int(settings, "mp_mode", 0);
    obs_data_set_default_double(settings, "mp_q", 0.70);
    obs_data_set_default_double(settings, "mp_wet_mix", 1.0);
    obs_data_set_default_double(settings, "mp_dry_mix", 0.0);
}

struct obs_source_info stir_midpass_info = {.id = "stir_midpass",
                                             .type = OBS_SOURCE_TYPE_FILTER,
                                             .output_flags = OBS_SOURCE_AUDIO,
                                             .get_name = stir_midpass_get_name,
                                             .create = stir_midpass_create,
                                             .destroy = stir_midpass_destroy,
                                             .filter_add = stir_midpass_add,
                                             .filter_remove = stir_midpass_remove,
                                             .get_properties = stir_midpass_properties,
                                             .get_defaults = stir_midpass_defaults,
                                             .update = stir_midpass_update};
