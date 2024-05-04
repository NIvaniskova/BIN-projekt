#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>

double get_threshold(std::vector<double> *desired, std::vector<double> *obtained) {
    std::vector<double>::iterator minimum = std::min_element((*obtained).begin(), (*obtained).end());

    double true_min = *minimum;

    std::vector<double>::iterator max = std::max_element((*obtained).begin(), (*obtained).end());

    double true_max = *max;

    double step = (true_max - true_min) / 10;

    //The true positive rate (TPR, also called sensitivity) is calculated as TP/TP+FN
    std::vector<double> thresholds;
    std::vector<double> tpr;
    std::vector<double> fpr;
    for (double threshold = true_min; threshold < true_max; threshold += step) {
        std::vector<int> class_obtained;
        for (int i = 0; i < (*obtained).size(); i++) {
            if ((*obtained)[i] < threshold) {
                class_obtained.push_back(0);
            } else {
                class_obtained.push_back(1);
            }
        }
        int tp = 0;
        int fn = 0;
        int fp = 0;
        int tn = 0;
        for (int j = 0; j < class_obtained.size(); j++) {
            if (class_obtained[j] == 1 and (*desired)[j] == 1) {
                tp++;
            } else if (class_obtained[j] == 0 and (*desired)[j] == 1) {
                fn++;
            } else if (class_obtained[j] == 1 and (*desired)[j] == 0) {
                fp++;
            } else {
                tn++;
            }
        }

        tpr.push_back((double) tp / (tp + fn));
        fpr.push_back((double) fp / (fp + tn));
        thresholds.push_back(threshold);
    }

    assert(!thresholds.empty());

    double best_threshold = thresholds[0];
    double max_tpr_fpr_diff = 0;
    for (int i = 0; i < tpr.size(); i++) {
        double diff = tpr[i] - fpr[i];
        if (diff > max_tpr_fpr_diff) {
            max_tpr_fpr_diff = diff;
            best_threshold = thresholds[i];
        }
    }
    return best_threshold;
}