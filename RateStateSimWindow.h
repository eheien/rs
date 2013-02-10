#include <qapplication.h>
#include <qmainwindow.h>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>

#ifndef _RATE_STATE_SIM_WINDOW_H_
#define _RATE_STATE_SIM_WINDOW_H_

class RateStateSimWindow : public QMainWindow {
	Q_OBJECT
	
private:
	double			param_a, param_b, param_k, param_r, param_w;
	QwtPlotCurve	*position_data, *velocity_data, *theta_data, *force_data, *driver_data;
	QwtPlot			*position_plot, *velocity_plot, *theta_plot, *force_plot;
	
public:
    RateStateSimWindow(QWidget *parent = 0, Qt::WindowFlags flags = 0);
	
public slots:
	void a_param_changed(const double &new_val) { param_a = pow(10,new_val); emit recalc(); };
	void ab_ratio_changed(const double &new_val) { param_b = param_a/new_val; emit recalc(); };
	void k_param_changed(const double &new_val) { param_k = new_val; emit recalc(); };
	void r_param_changed(const double &new_val) { param_r = pow(10,new_val); emit recalc(); };
	void recalc(void);
	
};

#endif
