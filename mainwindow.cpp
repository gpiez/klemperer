#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "displaywidget.h"

#include "integrator.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButtonStart_clicked()
{
    double centralMass = ui->centralMass->value();
    unsigned n = ui->nElements->value();

    integrator = new Integrator(n, centralMass);
    connect(integrator, SIGNAL(finished()), integrator, SLOT(deleteLater()));
    ui->displayWidget->setIntegrator(integrator);

    integrator->start();
    displayTimer = new QTimer(this);
    connect(displayTimer, SIGNAL(timeout()), this, SLOT(update()));
    displayTimer->start(16);
}

void MainWindow::on_pushButtonStop_clicked()
{
    if (integrator) {
        integrator->requestInterruption();
        integrator = 0;
        ui->displayWidget->setIntegrator(integrator);
        delete displayTimer;
    }
}

void MainWindow::update()
{
    ui->labelStep->setText(QString::number(integrator->getNSteps()));
    ui->labelTime->setText(QString::number(integrator->getDt_t(), 'f', 0));
    ui->labelEnergy->setText(QString::number(integrator->getEtot(), 'f', 6));
    ui->e_kin->setText(QString::number(integrator->getEkin(), 'f', 6));
    ui->e_pot->setText(QString::number(integrator->getEpot(), 'f', 6));
    ui->displayWidget->update();
}

Integrator *MainWindow::getIntegrator() const
{
    return integrator;
}

void MainWindow::on_lock_valueChanged(int arg1)
{
    ui->displayWidget->setLockBody(arg1);
}

void MainWindow::on_step_valueChanged(int arg1)
{
    integrator->setDt_param(pow(2.0, arg1));
}
