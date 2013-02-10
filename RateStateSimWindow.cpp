#include "RateStateSimWindow.h"
#include "RateState.h"

#include <iostream>

#include <qboxlayout.h>
#include <QHBoxLayout>
#include <qwt_slider.h>
#include <qwt_scale_engine.h>

RateStateSimWindow::RateStateSimWindow(QWidget *parent, Qt::WindowFlags flags) : QMainWindow(parent, flags) {
	QWidget				*main_widget;
	QVBoxLayout			*control_layout, *plot_layout;
	QHBoxLayout			*main_layout;
	int					i;
	QwtSlider			*a_param_slider, *ab_ratio_slider, *k_param_slider, *r_param_slider;
	
	param_a = 0.0625;
	param_b = 0.125;//param_a*2.84184246806358;//
	param_k = 20;
	param_r = 1e-5;
	param_w = 0.5;
	
	a_param_slider = new QwtSlider(0);
	a_param_slider->setScale(-2, 0);
	a_param_slider->setRange(-2, 0);
	a_param_slider->setValue(log10(param_a));
	a_param_slider->setScalePosition(QwtSlider::TopScale);
	
	ab_ratio_slider = new QwtSlider(0);
	ab_ratio_slider->setScale(0.001, 1.0);
	ab_ratio_slider->setRange(0.001, 1.0);
	ab_ratio_slider->setValue(param_a/param_b);
	ab_ratio_slider->setScalePosition(QwtSlider::TopScale);
	
	k_param_slider = new QwtSlider(0);
	k_param_slider->setScale(15, 25);
	k_param_slider->setRange(15, 25);
	k_param_slider->setValue(param_k);
	k_param_slider->setScalePosition(QwtSlider::TopScale);
	
	r_param_slider = new QwtSlider(0);
	r_param_slider->setScale(-7, 1);
	r_param_slider->setRange(-7, 1);
	r_param_slider->setValue(log10(param_r));
	r_param_slider->setScalePosition(QwtSlider::TopScale);
	
    connect(a_param_slider, SIGNAL(valueChanged(double)), this, SLOT(a_param_changed(double)));
    connect(ab_ratio_slider, SIGNAL(valueChanged(double)), this, SLOT(ab_ratio_changed(double)));
    connect(k_param_slider, SIGNAL(valueChanged(double)), this, SLOT(k_param_changed(double)));
    connect(r_param_slider, SIGNAL(valueChanged(double)), this, SLOT(r_param_changed(double)));
	
	control_layout = new QVBoxLayout;
	control_layout->addWidget(a_param_slider);
	control_layout->addWidget(ab_ratio_slider);
	control_layout->addWidget(k_param_slider);
	control_layout->addWidget(r_param_slider);
	
	position_data = new QwtPlotCurve(QwtText("X"));
	velocity_data = new QwtPlotCurve(QwtText("V"));
	theta_data = new QwtPlotCurve(QwtText("Theta"));
	force_data = new QwtPlotCurve(QwtText("Force"));
	driver_data = new QwtPlotCurve(QwtText("Driver Plate"));
	
	position_plot = new QwtPlot(QwtText("X"));
	velocity_plot = new QwtPlot(QwtText("Velocity"));
	velocity_plot->setAxisScaleEngine(QwtPlot::yLeft, new QwtLog10ScaleEngine);
	theta_plot = new QwtPlot(QwtText("Theta"));
	force_plot = new QwtPlot(QwtText("Force"));
	
	double min_t = 0, max_t = 160;
	position_plot->setAxisScale(QwtPlot::xBottom, min_t, max_t, 0);
	position_plot->setAxisScale(QwtPlot::yLeft, 100, 140, 0);
	velocity_plot->setAxisScale(QwtPlot::xBottom, min_t, max_t, 0);
	theta_plot->setAxisScale(QwtPlot::xBottom, min_t, max_t, 0);
	force_plot->setAxisScale(QwtPlot::xBottom, min_t, max_t, 0);
	
	position_data->attach(position_plot);
	driver_data->attach(position_plot);
	velocity_data->attach(velocity_plot);
	theta_data->attach(theta_plot);
	force_data->attach(force_plot);
	
	plot_layout = new QVBoxLayout;
	plot_layout->setContentsMargins(0, 0, 0, 0);
	
	plot_layout->addWidget(position_plot, 0);
	plot_layout->addWidget(velocity_plot, 0);
	plot_layout->addWidget(theta_plot, 0);
	plot_layout->addWidget(force_plot, 0);
	
	main_layout = new QHBoxLayout;
	main_layout->addLayout(control_layout);
	main_layout->addLayout(plot_layout);
	
	main_widget = new QWidget;
	main_widget->setLayout(main_layout);
	setCentralWidget(main_widget);
	
	emit recalc();
};

void RateStateSimWindow::recalc(void) {
	std::vector<std::vector<realtype> >	results;
	double					time_max = 200, time_step = 0.01;
	RSParams				params(NBLOCKS, NEQ, NPARAMS, time_step, time_max);
	unsigned int			i, npoints;
	double					xi, vi, hi;
	
	std::cerr << param_a << " " << param_b << " " << param_k << " " << param_r << " " << param_w << std::endl;
	for (i=0;i<params.num_blocks();++i) {
		params.param(i, A_PARAM) = RCONST(param_a);
		params.param(i, B_PARAM) = RCONST(param_b);
		params.param(i, K_PARAM) = RCONST(param_k);
		params.param(i, R_PARAM) = RCONST(param_r);
		params.param(i, W_PARAM) = RCONST(param_w);
		params.init_val(i, EQ_X) = RCONST(-10.0);
		params.init_val(i, EQ_V) = RCONST(1.0);
		params.init_val(i, EQ_H) = RCONST(1.0);
	}
	
	run_rate_state_sim(results, params);
	npoints = results.size();
	QVector<QPointF>		xdata(npoints), vdata(npoints), hdata(npoints), fdata(npoints), dp_data(npoints);
	double force_integral=0, old_t=0;
	
	for (i=0;i<results.size();++i) {
		xi = results[i][1];
		vi = results[i][2];
		hi = results[i][3];
		xdata[i] = QPointF(results[i][0], xi);
		vdata[i] = QPointF(results[i][0], vi);
		hdata[i] = QPointF(results[i][0], hi);
		fdata[i] = QPointF(results[i][0], F(0, vi, hi, params));
		dp_data[i] = QPointF(results[i][0], results[i][0]);
		force_integral += (results[i][0]-old_t)*F(0, vi, hi, params);
		old_t = results[i][0];
	}
	force_integral /= (results.size()/(time_step*time_max));
	
	//std::cerr << force_integral << std::endl;
	
	position_data->setData(new QwtPointSeriesData(xdata));
	velocity_data->setData(new QwtPointSeriesData(vdata));
	theta_data->setData(new QwtPointSeriesData(hdata));
	force_data->setData(new QwtPointSeriesData(fdata));
	driver_data->setData(new QwtPointSeriesData(dp_data));
	
	position_plot->replot();
	velocity_plot->replot();
	theta_plot->replot();
	force_plot->replot();
}
