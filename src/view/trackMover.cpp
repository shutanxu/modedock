/*
 * trackMover.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: lincong
 */

#include "trackMover.h"

trackMover::trackMover():
m_paused(false), m_pressed(false) {
	// TODO Auto-generated constructor stub
	m_translation = QVector2D();
	m_lastTime = QTime::currentTime();
}

void trackMover::push(const QPointF& p ) {
    m_translation = translation();
    m_pressed = true;
    m_lastTime = QTime::currentTime();
    m_lastPos = p;
    m_angularVelocity = 0.0f;
}

void trackMover::move(const QPointF& p ) {
//    printf("move !\n");
    if (!m_pressed)
        return;
    QTime currentTime = QTime::currentTime();
    int msecs = m_lastTime.msecsTo(currentTime);
    if (msecs <= 10)
        return;

    QLineF delta(m_lastPos, p);
   m_angularVelocity = delta.length() / msecs;
//    m_translation = QVector2D( p.x() - m_lastPos.x(),  p.y() - m_lastPos.y() );
    m_translation += QVector2D( p.x() - m_lastPos.x(),  p.y() - m_lastPos.y() );
//    m_translation *= 10.0;
    m_lastPos = p;
    m_lastTime = currentTime;
}

void trackMover::release(const QPointF& p )
{
    // Calling move() caused the rotation to stop if the framerate was too low.
    move(p);
    m_pressed = false;
}

void trackMover::start()
{
    m_lastTime = QTime::currentTime();
    m_paused = false;
}

void trackMover::stop()
{
    m_translation = translation();
    m_paused = true;
}

QVector2D trackMover::translation()  const {
	if (m_paused || m_pressed)
		return m_translation;
//	m_angularVelocity = 10.0f;
    QTime currentTime = QTime::currentTime();
    float angle = m_angularVelocity * m_lastTime.msecsTo(currentTime);
    QVector2D mt = m_translation + angle * m_translation;
    // printf("%7.3f%7.3f\n", mt.x(), mt.y() );
    return (mt);

//    float angle = m_angularVelocity * m_lastTime.msecsTo(currentTime);
//    return QQuaternion::fromAxisAndAngle(m_axis, angle) * m_rotation;
}

trackMover::~trackMover() {
	// TODO Auto-generated destructor stub
}

