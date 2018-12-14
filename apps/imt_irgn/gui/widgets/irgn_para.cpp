/****************************************************************************
**
** Copyright (C) 2011 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#include <QtGui>
#include <limits>

#include "irgn_para.h"

//! [0]
IRGN_para::IRGN_para()
{
  createIRGNSpinBoxes();
  createSimulBoxes();
  createButtons();

    QVBoxLayout *layout = new QVBoxLayout;
    spinBoxesGroup->setMinimumSize(250,310);
    simulationGroup->setMinimumSize(250,150);
    layout->addWidget(spinBoxesGroup);
    layout->addWidget(simulationGroup);
    layout->addLayout(buttonsLayout);

    setLayout(layout);

    //resize(400, 400);
    //this->setFixedSize(400,400);

    setWindowTitle(tr("IRGN Parameter"));

    //set init values
    _irgnpara.alpha0 = _sb_alpha0->value();
    _irgnpara.alpha_min = _sb_alpha_min->value();
    _irgnpara.alpha_q = _sb_alpha_q->value();
    _irgnpara.beta0 = _sb_beta0->value();
    _irgnpara.beta_min = _sb_beta_min->value();
    _irgnpara.beta_q = _sb_beta_q->value();
    _irgnpara.maxit = _sb_maxit->value();
    _irgnpara.tvits = _sb_tvits->value();
    _irgnpara.tvmax = _sb_tvmax->value();

    _irgn_simul.factor_x = _sb_factor_x->value();
    _irgn_simul.factor_y = _sb_factor_y->value();
    _irgn_simul.ref_lines = _sb_reflines->value();
    _irgn_simul.active = false;
}
//! [0]


//! [1]
void IRGN_para::createIRGNSpinBoxes()
{
    spinBoxesGroup = new QGroupBox(tr("IRGN Parameter"));

    QLabel* _lb_alpha0 = new QLabel(tr("alpha0:"));
    QLabel* _lb_alpha_min= new QLabel(tr("alpha min:"));
    QLabel* _lb_alpha_q= new QLabel(tr("alpha_q:"));
    QLabel* _lb_beta0= new QLabel(tr("beta0:"));
    QLabel* _lb_beta_min= new QLabel(tr("beta min:"));
    QLabel* _lb_beta_q= new QLabel(tr("beta_q:"));
    QLabel* _lb_tvmax= new QLabel(tr("tv max:"));
    QLabel* _lb_tvits= new QLabel(tr("tv its:"));
    QLabel* _lb_maxit= new QLabel(tr("max it:"));

    _sb_alpha0 = new QDoubleSpinBox;    // initial penalty alpha_0 (L2, sensitivites)
    _sb_alpha_min = new QDoubleSpinBox; // final value of alpha
    _sb_alpha_q = new QDoubleSpinBox;   // reduction factor for alpha
    _sb_beta0 = new QDoubleSpinBox;     // initial penalty beta_0 (image)
    _sb_beta_min = new QDoubleSpinBox;  // final value of beta: 0: no T(G)V effect, >0 effect
    _sb_beta_q = new QDoubleSpinBox;    // reduction factor for beta

    _sb_tvmax = new QSpinBox;     // upper bound on number of gradient steps
    _sb_tvits = new QSpinBox;     // initial number of gradient steps
    _sb_maxit = new QSpinBox;     // maximum number of IRGN iterations


    _sb_tvmax->setRange(0, 10000);
    _sb_tvmax->setSingleStep(10);
    _sb_tvmax->setValue(1000);

    _sb_tvits->setRange(0, 1000);
    _sb_tvits->setSingleStep(1);
    _sb_tvits->setValue(20);

    _sb_maxit->setRange(0, 100);
    _sb_maxit->setSingleStep(1);
    _sb_maxit->setValue(5);

    _sb_alpha0->setRange(std::numeric_limits<float>::min(), 100);  //1.17549e-38
    _sb_alpha0->setSingleStep(0.005);
    _sb_alpha0->setValue(1);
    _sb_alpha0->setDecimals(3);

    _sb_alpha_min->setRange(std::numeric_limits<float>::min(), 100);
    _sb_alpha_min->setSingleStep(0.005);
    _sb_alpha_min->setValue(std::numeric_limits<float>::min());
    _sb_alpha_min->setDecimals(3);

    _sb_alpha_q->setRange(std::numeric_limits<float>::min(), 100);
    _sb_alpha_q->setSingleStep(0.005);
    _sb_alpha_q->setValue(double(0.1));
    _sb_alpha_q->setDecimals(3);

    _sb_beta0->setRange(std::numeric_limits<float>::min(), 100);
    _sb_beta0->setSingleStep(0.005);
    _sb_beta0->setValue(1);
    _sb_beta0->setDecimals(3);

    _sb_beta_min->setRange(std::numeric_limits<float>::min(), 100);
    _sb_beta_min->setSingleStep(0.005);
    _sb_beta_min->setValue(std::numeric_limits<float>::min());
    _sb_beta_min->setDecimals(3);

    _sb_beta_q->setRange(std::numeric_limits<float>::min(), 100);
    _sb_beta_q->setSingleStep(0.005);
    _sb_beta_q->setValue(double(0.2));
    _sb_beta_q->setDecimals(3);


    QGridLayout *spinBoxLayout = new QGridLayout;
    spinBoxLayout->addWidget(_lb_alpha0, 0, 0, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_lb_alpha_min, 1, 0, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_lb_alpha_q, 2, 0, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_lb_beta0, 3, 0, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_lb_beta_min, 4, 0, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_lb_beta_q, 5, 0, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_lb_tvmax, 6, 0, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_lb_tvits, 7, 0, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_lb_maxit, 8, 0, 2, 2, Qt::AlignLeft);

    spinBoxLayout->addWidget(_sb_alpha0, 0, 1, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_sb_alpha_min, 1, 1, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_sb_alpha_q, 2, 1, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_sb_beta0, 3, 1, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_sb_beta_min, 4, 1, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_sb_beta_q, 5, 1, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_sb_tvmax, 6, 1, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_sb_tvits, 7, 1, 2, 2, Qt::AlignLeft);
    spinBoxLayout->addWidget(_sb_maxit, 8, 1, 2, 2, Qt::AlignLeft);

    spinBoxesGroup->setLayout(spinBoxLayout);
}

//! [1]
void IRGN_para::createSimulBoxes()
{
    simulationGroup = new QGroupBox(tr("Simulation Parameter"));

    QLabel* _lb_factor_x = new QLabel(tr("factor_x:"));
    QLabel* _lb_factor_y= new QLabel(tr("factor_y:"));
    QLabel* _lb_ref_lines= new QLabel(tr("ref-lines:"));

    _sb_factor_x = new QSpinBox;    // simul factor_x
    _sb_factor_y = new QSpinBox;    // simul factor_y
    _sb_reflines = new QSpinBox;    // simul reflines

    _simul_active = new QCheckBox;



    _sb_factor_x->setRange(1, 50);
    _sb_factor_x->setSingleStep(1);
    _sb_factor_x->setValue(1);

    _sb_factor_y->setRange(1, 50);
    _sb_factor_y->setSingleStep(1);
    _sb_factor_y->setValue(1);

    _sb_reflines->setRange(0, 200);
    _sb_reflines->setSingleStep(1);
    _sb_reflines->setValue(0);



    QGridLayout *simuLayout = new QGridLayout;
    simuLayout->addWidget(_lb_factor_x, 0, 0, 2, 2, Qt::AlignLeft);
    simuLayout->addWidget(_lb_factor_y, 1, 0, 2, 2, Qt::AlignLeft);
    simuLayout->addWidget(_lb_ref_lines, 2, 0, 2, 2, Qt::AlignLeft);

    simuLayout->addWidget(_sb_factor_x, 0, 1, 2, 2, Qt::AlignLeft);
    simuLayout->addWidget(_sb_factor_y, 1, 1, 2, 2, Qt::AlignLeft);
    simuLayout->addWidget(_sb_reflines, 2, 1, 2, 2, Qt::AlignLeft);

    simulationGroup->setLayout(simuLayout);
}


//! [1]
void IRGN_para::createButtons()
{
  QPushButton *pb_ok = createButton("Apply", SLOT(okbutton()));
  QPushButton *pb_cancle = createButton("Close", SLOT(cancelbutton()));

  buttonsLayout = new QHBoxLayout;
  //buttonsLayout->addStretch();
  buttonsLayout->addWidget(pb_ok);
  buttonsLayout->addWidget(pb_cancle);

}


QPushButton *IRGN_para::createButton(const QString &text, const char *member)
{
    QPushButton *button = new QPushButton(text);
    connect(button, SIGNAL(clicked()), this, member);
    return button;
}

void IRGN_para::cancelbutton()
{
  _sb_alpha0->setValue(_irgnpara.alpha0);
  _sb_alpha_min->setValue(_irgnpara.alpha_min);
  _sb_alpha_q->setValue(_irgnpara.alpha_q);
  _sb_beta0->setValue(_irgnpara.beta0);
  _sb_beta_min->setValue(_irgnpara.beta_min);
  _sb_beta_q->setValue(_irgnpara.beta_q);
  _sb_maxit->setValue(_irgnpara.maxit);
  _sb_tvits->setValue(_irgnpara.tvits);
  _sb_tvmax->setValue(_irgnpara.tvmax);

  _sb_factor_x->setValue(_irgn_simul.factor_x);
  _sb_factor_y->setValue(_irgn_simul.factor_y);
  _sb_reflines->setValue(_irgn_simul.ref_lines);

  this->close();
}


void IRGN_para::okbutton()
{
  _irgnpara.alpha0 = _sb_alpha0->value();
  _irgnpara.alpha_min = _sb_alpha_min->value();
  _irgnpara.alpha_q = _sb_alpha_q->value();
  _irgnpara.beta0 = _sb_beta0->value();
  _irgnpara.beta_min = _sb_beta_min->value();
  _irgnpara.beta_q = _sb_beta_q->value();
  _irgnpara.maxit = _sb_maxit->value();
  _irgnpara.tvits = _sb_tvits->value();
  _irgnpara.tvmax = _sb_tvmax->value();

  _irgn_simul.factor_x = _sb_factor_x->value();
  _irgn_simul.factor_y = _sb_factor_y->value();
  _irgn_simul.ref_lines = _sb_reflines->value();
  if((_irgn_simul.factor_x > 1) || (_irgn_simul.factor_y > 1))
    _irgn_simul.active = true;
  else
    _irgn_simul.active = false;

  //this->close();
}
