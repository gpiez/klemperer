#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>

class Integrator;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    Integrator *getIntegrator() const;

private slots:
    void on_pushButtonStart_clicked();

    void on_pushButtonStop_clicked();
    void update();
    void on_lock_valueChanged(int arg1);

    void on_step_valueChanged(int arg1);

private:
    Ui::MainWindow *ui;
    Integrator* integrator;
    QTimer* displayTimer;
};

#endif // MAINWINDOW_H
