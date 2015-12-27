#ifndef DISPLAYWIDGET_H
#define DISPLAYWIDGET_H

#include <QMatrix4x4>
//#include <QVector2D>
#include <QVector3D>
#include <QOpenGLBuffer>
#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions>

#include <QOpenGLDebugLogger>

class Integrator;

struct VertexData
{
    QVector3D position;
};

class DisplayWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT
public:
    explicit DisplayWidget(QWidget *parent);
    ~DisplayWidget();

    void setIntegrator(Integrator *value);

    void setLockBody(int value);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void wheelEvent(QWheelEvent* ) override;

private:
    void initShaders();
    void initTextures();
    void check();

    double scale;

    int modelViewLocation;
    int sphereRadiusLocation;
    int sphereColorLocation;
    int lightPositionLocation;
    int positionLocation;
    int lockBody;

    QVector3D lightPosition;
    QVector3D color[101];
    QOpenGLBuffer array;
    QOpenGLBuffer index;
    QOpenGLShaderProgram sphere;
    QOpenGLDebugLogger* logger;

    Integrator* integrator;
};

#endif // DISPLAYWIDGET_H
