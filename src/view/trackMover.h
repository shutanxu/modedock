/*
 * trackMover.h
 *
 *  Created on: Sep 4, 2014
 *      Author: lincong
 */

#ifndef TRACKMOVER_H_
#define TRACKMOVER_H_
#include <QtGui>
#include <QtGui/QVector2D>

class trackMover {
public:
	trackMover();
	virtual ~trackMover();

	void push(const QPointF& p);
	void move(const QPointF& p );
	void release(const QPointF& p);
	void start(); // starts clock
	void stop(); // stops clock
	QVector2D translation() const ;
	QPointF getMousePrevPos () const { return m_lastPos;};

	private:
		float m_angularVelocity;
		QVector2D m_translation;
    	// QVector3D m_axis;
	    QPointF m_lastPos;
	    QTime m_lastTime;
	    bool m_paused;
	    bool m_pressed;
};

#endif /* TRACKMOVER_H_ */
