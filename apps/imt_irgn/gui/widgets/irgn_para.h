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

#ifndef IRGN_para_H
#define IRGN_para_H

#include <QWidget>
//#include "agile/calc/irgn_all.hpp"
#include "agile/calc/irgn.hpp"

QT_BEGIN_NAMESPACE
class QDateTimeEdit;
class QDoubleSpinBox;
class QSpinBox;
class QGroupBox;
class QLabel;
class QHBoxLayout;
class QPushButton;
class QCheckBox;
QT_END_NAMESPACE

struct IRGN_Simul{
  bool active;
  unsigned factor_x;
  unsigned factor_y;
  unsigned ref_lines;
};


//! [0]
class IRGN_para : public QWidget
{
    Q_OBJECT

public:
    IRGN_para();

    agile::IRGN_Params getIRGN_params()
    {
      return _irgnpara;
    }

    IRGN_Simul getIRGN_simul()
    {
      return _irgn_simul;
    }

public slots:
    void okbutton();
    void cancelbutton();

private:
    QPushButton *createButton(const QString &text, const char *member);

    void createIRGNSpinBoxes();
    void createSimulBoxes();
    void createButtons();

    void createSpinBoxes();
    void createDoubleSpinBoxes();

    QGroupBox *spinBoxesGroup;
    QGroupBox *simulationGroup;
    QHBoxLayout *buttonsLayout;

    QDoubleSpinBox* _sb_alpha0;        // initial penalty alpha_0 (L2, sensitivites)
    QDoubleSpinBox* _sb_alpha_min;     // final value of alpha
    QDoubleSpinBox* _sb_alpha_q;       // reduction factor for alpha
    QDoubleSpinBox* _sb_beta0;         // initial penalty beta_0 (image)
    QDoubleSpinBox* _sb_beta_min;      // final value of beta: 0: no T(G)V effect, >0 effect
    QDoubleSpinBox* _sb_beta_q;        // reduction factor for beta

    QSpinBox* _sb_factor_x;        // simul factor_x
    QSpinBox* _sb_factor_y;        // simul factor_y
    QSpinBox* _sb_reflines;        // simul reflines
    QCheckBox*      _simul_active;       // simul active

    QSpinBox *_sb_tvmax;     // upper bound on number of gradient steps
    QSpinBox *_sb_tvits;     // initial number of gradient steps
    QSpinBox *_sb_maxit;     // maximum number of IRGN iterations



    agile::IRGN_Params _irgnpara;
    IRGN_Simul _irgn_simul;

};
//! [0]

#endif
