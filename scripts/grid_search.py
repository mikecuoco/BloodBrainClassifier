from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

# model = sklearn.svm.SVC()
# tuned_parameters = [{'kernel': ['poly'], 'gamma': [1e-7, 1e-8, 1e-10], 'C': [1, 0.000001], 'degree':[5, 6, 7, 8, 9]},
#                     {'kernel': ['rbf'], 'gamma': [1e-3], 'C': [1]},
#                     {'kernel': ['linear'], 'C': [1]}]

def main(data, model, tuned_parameters, score = 'f1', k=4):
    """
        input:
            data - a tuple of X_train, y_train, X_test, y_test
            model = an sklearn classifier, like SVC()
            tuned_parameters - a list of dictionaries, where each dictionary maps a hyperparameter to a list of possible values
        output:
            clf - the GridSearch classifier object
    """
    X_train, y_train, X_test, y_test = data
    print("# Tuning hyper-parameters for %s" % score)
    print()

    clf = GridSearchCV(
        model, tuned_parameters, scoring='%s_macro' % score, cv=k
    )
    clf.fit(X_train, y_train)

    print("Best parameters set found on development set:")
    print()
    print(clf.best_params_)
    print()
    print("Grid scores on development set:")
    print()
    means = clf.cv_results_['mean_test_score']
    stds = clf.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))
    print()

    print("Detailed classification report:")
    print()
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print()
    y_true, y_pred = y_test, clf.predict(X_test)
    print(classification_report(y_true, y_pred))
    print()

    return clf
