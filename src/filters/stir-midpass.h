#pragma once

#include "obs.h"

const char *stir_midpass_get_name(void *data);
void *stir_midpass_create(obs_data_t *settings, obs_source_t *source);
extern struct obs_source_info stir_midpass_info;
