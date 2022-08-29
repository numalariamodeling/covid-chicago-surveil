#ifndef TIME_SERIES_H
#define TIME_SERIES_H
#include "Utility.h"

struct TimeSeriesAnchorPoint {
    size_t sim_day;
    double value;
};

vector<double> linInterpolate(double start, double end, size_t n) {
    // "end" is the value immediately after the sequence produced here
    // e.g. linInterpolate(0.0, 1.0, 5) produces {0.0, 0.2, 0.4, 0.6, 0.8}
    assert(n>=1);
    //cout << "Start: " << start << "; End: " << end << endl;
    
    // initialize the n values to the start value
    double step = (end - start)/n;
    vector<double> v(n, start);
    for (size_t i = 1; i < n; i++) { v[i] += i * step; }
    return v;
}

vector<double> linInterpolateTimeSeries(vector<TimeSeriesAnchorPoint> ap) {
    if (ap.size() < 2) {
        cerr << "Must have at least two anchor points." << endl;
        exit(1);
    }
    sort(ap.begin(), ap.end(), [](TimeSeriesAnchorPoint &a, TimeSeriesAnchorPoint &b){return a.sim_day < b.sim_day;});
    vector<double> v;

    for (size_t i = 0; i < ap.size() - 1; i++) {
        const TimeSeriesAnchorPoint A = ap[i];
        const TimeSeriesAnchorPoint B = ap[i+1];
        size_t n = B.sim_day - A.sim_day;

        vector<double> v_i = linInterpolate(A.value, B.value, n);
        v.insert(v.end(), v_i.begin(), v_i.end());
    }
    v.push_back(ap.back().value); // bc linInterpolate leaves off last value, to prevent repetitions

    //int series_julian_start = to_sim_day(julian_start_year, julian_start_day, ap[0].date);
    //if (series_julian_start < 0) {
    //    v.erase(v.begin(), v.begin()-series_julian_start);
    //}

    return v;
}

vector<double> stepwiseTimeSeries(vector<TimeSeriesAnchorPoint> ap) {
    if (ap.size() < 2) {
        cerr << "Must have at least two anchor points." << endl;
        exit(1);
    }
    sort(ap.begin(), ap.end(), [](TimeSeriesAnchorPoint &a, TimeSeriesAnchorPoint &b){return a.sim_day < b.sim_day;});
    vector<double> v;

    for (size_t i = 0; i < ap.size() - 1; i++) {
        const TimeSeriesAnchorPoint A = ap[i];
        const TimeSeriesAnchorPoint B = ap[i+1];
        size_t n = B.sim_day - A.sim_day;
        
        for (size_t j = 0; j < n; j++) {
            v.push_back(A.value);
        }
    }
    v.push_back(ap.back().value); // bc linInterpolate leaves off last value, to prevent repetitions
    return v;
}
#endif
