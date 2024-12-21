#include <ooo_cpu.h>
#include "tage.h"

void O3_CPU::initialize_branch_predictor()
{
    tage_init();
}

uint8_t O3_CPU::predict_branch(uint64_t ip)
{
    return tage_predict(ip);
}

void O3_CPU::last_branch_result(uint64_t ip, uint64_t branch_target, uint8_t taken, uint8_t branch_type)
{
    tage_train(ip,taken);
}