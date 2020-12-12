#include <paradem/timeInfo.h>

TimeInfo::TimeInfo() {
    calc = overall = io = 0;
    vmpeak = vmhwm = 0;
}

TimeInfo::TimeInfo(double calc, double overall, double io, long vmpeak, long vmhwml) : calc(calc), overall(overall), io(io), vmpeak(vmpeak), vmhwm(vmhwml) {}
