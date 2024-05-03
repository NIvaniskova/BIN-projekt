#include "cgp.h"
#include "symb.h"
#include <vector>
#include <algorithm>

double get_threshold( std::vector<double> desired,  std::vector<double> obtained){
    auto min = std::min_element(obtained.begin(), obtained.end());

    double true_min;
    if (min == obtained.end()) {
        return -1; // List was empty //TODO
    } else {
        true_min = min - obtained.begin();
    }

    auto max = std::max_element(obtained.begin(), obtained.end());

    double true_max;
    if (max == obtained.end()) {
        return -1; // List was empty //TODO
    } else {
        true_max = max - obtained.begin();
    }

    double step = (true_max - true_min) / 10;

    //The true positive rate (TPR, also called sensitivity) is calculated as TP/TP+FN

    std::vector<double> thresholds;
    std::vector<double> tpr;
    std::vector<double> fpr;
    for (double threshold = true_min; threshold < true_max; threshold + step){
        std::vector<int> class_obtained;
        for (int i = 0; i < obtained.size(); i++){
            if (obtained[i] < threshold){
                class_obtained[i] = 0;
            }else{
                class_obtained[i] = 1;
            }
        }
        int tp;
        int fn;
        int fp;
        int tn;
        for (int j = 0; j < class_obtained.size(); j++){
            if (class_obtained[j] == 1 and desired[j] == 1){
                tp++;
            }else if (class_obtained[j] == 0 and desired[j] == 1){
                fn++;
            }else if (class_obtained[j] == 1 and desired[j] == 0){
                fp++;
            }else{
                tn++;
            }
        }
        tpr.push_back(tp/(tp+fn));
        fpr.push_back(fp/(fp+tn));
        thresholds.push_back(threshold);
    }

    double best_threshold = thresholds[0];
    double max_tpr_fpr_diff = 0;
    for(int i = 0; i < tpr.size(); i++){
        double diff = tpr[i] - fpr[i];
        if ( diff > max_tpr_fpr_diff){
            max_tpr_fpr_diff = diff;
            best_threshold = thresholds[i];
        }
    }
    return best_threshold;
}