#include "displaywidget.h"

#include <QOpenGLFunctions>
#include <QWheelEvent>
#include <QtDebug>
#include <integrator.h>

#define CHECK do {\
    QList<QOpenGLDebugMessage> messages = logger->loggedMessages();\
    foreach (const QOpenGLDebugMessage &message, messages) {\
        qDebug() << __FILE__ << __LINE__ << message.message();\
   } } while (0);

#define CHECK2(x) if ((x)<0) { qDebug() << __FILE__ << __LINE__ << " invalid value"; return; }

DisplayWidget::DisplayWidget(QWidget *parent) :
    QOpenGLWidget(parent),
    scale(1.0),
    index(QOpenGLBuffer::IndexBuffer),
    logger(new QOpenGLDebugLogger(this)),
    integrator(0)
{ }

DisplayWidget::~DisplayWidget()
{
    makeCurrent();
    delete logger;
    doneCurrent();
}

void DisplayWidget::initializeGL()
{
    initializeOpenGLFunctions();
    logger->initialize();
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    initShaders();
    initTextures();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    VertexData vertices[] = {
        { QVector3D(0, 0, 0) },
        { QVector3D(0, 0, 0) },
        { QVector3D(0, 0, 0) },
        { QVector3D(0, 0, 0) }
    };
    array.create(); CHECK;
    array.bind(); CHECK;
    array.allocate(vertices, sizeof(vertices)); CHECK;

    GLushort indices[] = {
           0,  1,  2,  3
    };
    index.create(); CHECK;
    index.bind(); CHECK;
    index.allocate(indices, sizeof(indices)); CHECK;

    modelViewLocation = sphere.uniformLocation("modelView"); CHECK;
    sphereRadiusLocation = sphere.uniformLocation("sphereRadius"); CHECK;
    sphereColorLocation = sphere.uniformLocation("sphereColor"); CHECK;
    lightPositionLocation = sphere.uniformLocation("lightPosition"); CHECK;
    positionLocation = sphere.attributeLocation("position"); CHECK;

    lightPosition = QVector3D(5,2,10);
    lightPosition.normalize();

    for (int i=1; i<=100; i++) {
        color[i] = QVector3D(fmod(0.8321f*i,1), fmod(0.5195f*i,1), fmod(0.7123f*i,1));
        color[i].normalize();
    }
    color[0]=QVector3D(1.0,1.0,1.0);
}

void DisplayWidget::resizeGL(int w, int h)
{
    // Update projection matrix and other size related settings:
    QMatrix4x4 projection;
    projection.setToIdentity();
    projection.perspective(90.0f, w / float(h), 1.0f, 10000.0f);
    sphere.setUniformValue("projection", projection); CHECK;
}

void DisplayWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT);
    if (integrator) {
        QMatrix4x4 rotateView;
        if (lockBody > 0) {
            double d1 = integrator->getPos(lockBody, 0)-integrator->getPos(0, 0);
            double d2 = integrator->getPos(lockBody, 1)-integrator->getPos(0, 1);
            QVector3D axis = QVector3D(0, 0, -1);
            float angle = 180.0/M_PI*atan2(d2,d1);
            rotateView.rotate(angle, axis);
        }
        for (unsigned i=0; i<integrator->getN(); i++) {
            QMatrix4x4 modelView;
            modelView *= rotateView;
            modelView.translate(integrator->getPos(i, 0),
                                integrator->getPos(i, 1),
                                integrator->getPos(i, 2)-400*scale);
            sphere.setUniformValue(modelViewLocation, modelView); CHECK;

            sphere.setUniformValue(sphereRadiusLocation, float(integrator->getRadius(i))); CHECK;
            sphere.setUniformValue(sphereColorLocation, color[i]); CHECK;

            sphere.setUniformValue(lightPositionLocation, lightPosition); CHECK;
            // Draw the scene:
            array.bind(); CHECK;
            index.bind(); CHECK;

            int vertexPosition = sphere.attributeLocation("position"); CHECK; CHECK2(vertexPosition);
            sphere.enableAttributeArray(vertexPosition); CHECK;
            sphere.setAttributeBuffer(vertexPosition, GL_FLOAT, 0, 3); CHECK;

            glDrawElements(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_SHORT, 0); CHECK;
        }
    }

}

void DisplayWidget::wheelEvent(QWheelEvent* event)
{
    scale *= (1.0 + event->delta()/1200.0);
}

void DisplayWidget::initShaders()
{
    if (!sphere.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/sphere.vsh")) {
        qCritical() << sphere.log();
        return;
    }
    if (!sphere.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/sphere.fsh")) {
        qCritical() << sphere.log();
        return;
    }
    if (!sphere.link()) {
        qCritical() << sphere.log();
        return;
    }
    if (!sphere.bind()) {
        qCritical() << sphere.log();
        return;
    }
    qDebug() << sphere.log();

}

void DisplayWidget::initTextures()
{

}

void DisplayWidget::setLockBody(int value)
{
    lockBody = value;
}

void DisplayWidget::setIntegrator(Integrator *value)
{
    integrator = value;
}

