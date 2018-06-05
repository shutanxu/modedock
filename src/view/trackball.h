/**
 *
If you are using Orthographic projection then your object looks the same size irrespective of
 your camera position. So here if you want to zoom you should modify the viewing volume
  by adjusting the parameters to the glOrtho. Increasing the volume will make the object
   look smaller and decreasing the volume will make the object look bigger. This is more a
   natural way. I have seen some people use the glScale too but the previous one is more
   natural in my opinion.

If you want to do zooming in perspective projection then @David M suggestion of modifying
 the FOV is fine.
 *
 */
#ifndef TRACKBALL_H
#define TRACKBALL_H

#include <QtGui>

#include <QtGui/qvector3d.h>
#include <QtGui/qquaternion.h>
#include<iostream>
using namespace std;
class TrackBall
{
public:
    enum TrackMode
    {
        Plane,
        Sphere,
    };
    TrackBall(TrackMode mode = Sphere);
    TrackBall(float angularVelocity, const QVector3D& axis, TrackMode mode = Sphere);
    // coordinates in [-1,1]x[-1,1]
    void push(const QPointF& p, const QQuaternion &transformation);
    void move(const QPointF& p, const QQuaternion &transformation);
    void release(const QPointF& p, const QQuaternion &transformation);
    void start(); // starts clock
    void stop(); // stops clock
    QQuaternion rotation() const;
    QPointF getMousePrevPos () const { return m_lastPos;};

private:
    QQuaternion m_rotation;
    QVector3D m_axis;
    float m_angularVelocity;

    QPointF m_lastPos;
    QTime m_lastTime;
    bool m_paused;
    bool m_pressed;
    TrackMode m_mode;
};

#endif
