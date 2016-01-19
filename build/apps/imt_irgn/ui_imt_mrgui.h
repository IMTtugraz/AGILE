/********************************************************************************
** Form generated from reading UI file 'imt_mrgui.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_IMT_MRGUI_H
#define UI_IMT_MRGUI_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QFormLayout>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QStatusBar>
#include <QtGui/QTextBrowser>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_IMT_MRGui
{
public:
    QAction *actionRead_DicomRaw;
    QAction *actionOpen_meas_dat;
    QAction *actionParameter;
    QAction *actionOpen_matlab_bin;
    QAction *actionAuto_Archive;
    QAction *actionSpecial;
    QWidget *centralWidget;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_2;
    QFormLayout *formLayout_stdinfo;
    QLabel *label_rows;
    QLabel *label_columns;
    QLabel *label_coils;
    QLabel *label_rowsval;
    QLabel *label_columnsval;
    QLabel *label_coilsval;
    QFrame *frame_3;
    QFrame *frame_4;
    QWidget *gridLayoutWidget;
    QGridLayout *gridLayout_imginfo;
    QLabel *label_slival;
    QLabel *label_slices;
    QLabel *label_acq;
    QLabel *label_acqval;
    QLabel *label_rep;
    QLabel *label_repval;
    QLabel *label_par;
    QLabel *label_parval;
    QLabel *label_echo;
    QLabel *label_echoval;
    QLabel *label_pha;
    QLabel *label_phaval;
    QLabel *label_set;
    QLabel *label_setval;
    QLabel *label_seg;
    QLabel *label_segval;
    QGroupBox *groupBox;
    QGridLayout *gridLayout;
    QVBoxLayout *verticalLayout;
    QVBoxLayout *verticalLayout_3;
    QCheckBox *cb_cropFOV;
    QRadioButton *radioButton_rawdata;
    QRadioButton *radioButton_ifft;
    QRadioButton *radioButton_irgnl2;
    QRadioButton *radioButton_irgntv;
    QRadioButton *radioButton_irgntgv;
    QFrame *frame;
    QGridLayout *gridLayout_imgselect;
    QLabel *label_slices_2;
    QLabel *label_acq_2;
    QLabel *label_rep_2;
    QLabel *label_echo_2;
    QLabel *label_pha_2;
    QLabel *label_set_2;
    QLabel *label_seg_2;
    QComboBox *combo_sli_select;
    QComboBox *combo_acq_select;
    QComboBox *combo_rep_select;
    QComboBox *combo_par_select;
    QComboBox *combo_echo_select;
    QComboBox *combo_pha_select;
    QComboBox *combo_seg_select;
    QComboBox *combo_set_select;
    QLabel *label_par_2;
    QFrame *frame_2;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout;
    QGroupBox *groupBox_3;
    QWidget *verticalLayoutWidget_2;
    QVBoxLayout *verticalLayout_2;
    QTextBrowser *text_info;
    QVBoxLayout *verticalLayout_4;
    QSpacerItem *verticalSpacer;
    QGroupBox *groupBox_5;
    QWidget *formLayoutWidget;
    QFormLayout *formLayout;
    QPushButton *pb_outputpath;
    QPushButton *pb_automatic;
    QLabel *label;
    QLabel *label_2;
    QComboBox *cb_savingtype;
    QLabel *label_3;
    QGroupBox *groupBox_4;
    QWidget *formLayoutWidget_3;
    QFormLayout *formLayout_2;
    QPushButton *pushButton;
    QCheckBox *cb_autofilename;
    QProgressBar *progressBar;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuIRGN_Config;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *IMT_MRGui)
    {
        if (IMT_MRGui->objectName().isEmpty())
            IMT_MRGui->setObjectName(QString::fromUtf8("IMT_MRGui"));
        IMT_MRGui->resize(687, 571);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(IMT_MRGui->sizePolicy().hasHeightForWidth());
        IMT_MRGui->setSizePolicy(sizePolicy);
        actionRead_DicomRaw = new QAction(IMT_MRGui);
        actionRead_DicomRaw->setObjectName(QString::fromUtf8("actionRead_DicomRaw"));
        actionOpen_meas_dat = new QAction(IMT_MRGui);
        actionOpen_meas_dat->setObjectName(QString::fromUtf8("actionOpen_meas_dat"));
        actionParameter = new QAction(IMT_MRGui);
        actionParameter->setObjectName(QString::fromUtf8("actionParameter"));
        actionOpen_matlab_bin = new QAction(IMT_MRGui);
        actionOpen_matlab_bin->setObjectName(QString::fromUtf8("actionOpen_matlab_bin"));
        actionAuto_Archive = new QAction(IMT_MRGui);
        actionAuto_Archive->setObjectName(QString::fromUtf8("actionAuto_Archive"));
        actionSpecial = new QAction(IMT_MRGui);
        actionSpecial->setObjectName(QString::fromUtf8("actionSpecial"));
        centralWidget = new QWidget(IMT_MRGui);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        groupBox_2 = new QGroupBox(centralWidget);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        groupBox_2->setGeometry(QRect(10, 10, 631, 101));
        groupBox_2->setAutoFillBackground(false);
        gridLayout_2 = new QGridLayout(groupBox_2);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setContentsMargins(11, 11, 11, 11);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        formLayout_stdinfo = new QFormLayout();
        formLayout_stdinfo->setSpacing(6);
        formLayout_stdinfo->setObjectName(QString::fromUtf8("formLayout_stdinfo"));
        formLayout_stdinfo->setSizeConstraint(QLayout::SetFixedSize);
        label_rows = new QLabel(groupBox_2);
        label_rows->setObjectName(QString::fromUtf8("label_rows"));

        formLayout_stdinfo->setWidget(0, QFormLayout::LabelRole, label_rows);

        label_columns = new QLabel(groupBox_2);
        label_columns->setObjectName(QString::fromUtf8("label_columns"));

        formLayout_stdinfo->setWidget(1, QFormLayout::LabelRole, label_columns);

        label_coils = new QLabel(groupBox_2);
        label_coils->setObjectName(QString::fromUtf8("label_coils"));

        formLayout_stdinfo->setWidget(2, QFormLayout::LabelRole, label_coils);

        label_rowsval = new QLabel(groupBox_2);
        label_rowsval->setObjectName(QString::fromUtf8("label_rowsval"));
        label_rowsval->setMaximumSize(QSize(70, 16777215));

        formLayout_stdinfo->setWidget(0, QFormLayout::FieldRole, label_rowsval);

        label_columnsval = new QLabel(groupBox_2);
        label_columnsval->setObjectName(QString::fromUtf8("label_columnsval"));
        label_columnsval->setMinimumSize(QSize(70, 0));
        label_columnsval->setMaximumSize(QSize(70, 16777215));

        formLayout_stdinfo->setWidget(1, QFormLayout::FieldRole, label_columnsval);

        label_coilsval = new QLabel(groupBox_2);
        label_coilsval->setObjectName(QString::fromUtf8("label_coilsval"));
        label_coilsval->setMinimumSize(QSize(70, 0));
        label_coilsval->setMaximumSize(QSize(7, 16777215));

        formLayout_stdinfo->setWidget(2, QFormLayout::FieldRole, label_coilsval);


        gridLayout_2->addLayout(formLayout_stdinfo, 0, 0, 2, 1);

        frame_3 = new QFrame(groupBox_2);
        frame_3->setObjectName(QString::fromUtf8("frame_3"));
        frame_3->setMaximumSize(QSize(10, 16777215));
        frame_3->setFrameShape(QFrame::NoFrame);
        frame_3->setFrameShadow(QFrame::Raised);

        gridLayout_2->addWidget(frame_3, 0, 2, 2, 1);

        frame_4 = new QFrame(groupBox_2);
        frame_4->setObjectName(QString::fromUtf8("frame_4"));
        frame_4->setMinimumSize(QSize(436, 0));
        frame_4->setFrameShape(QFrame::NoFrame);
        frame_4->setFrameShadow(QFrame::Raised);
        gridLayoutWidget = new QWidget(frame_4);
        gridLayoutWidget->setObjectName(QString::fromUtf8("gridLayoutWidget"));
        gridLayoutWidget->setGeometry(QRect(0, 0, 431, 61));
        gridLayout_imginfo = new QGridLayout(gridLayoutWidget);
        gridLayout_imginfo->setSpacing(0);
        gridLayout_imginfo->setContentsMargins(11, 11, 11, 11);
        gridLayout_imginfo->setObjectName(QString::fromUtf8("gridLayout_imginfo"));
        gridLayout_imginfo->setSizeConstraint(QLayout::SetMinimumSize);
        gridLayout_imginfo->setContentsMargins(0, 0, 0, 0);
        label_slival = new QLabel(gridLayoutWidget);
        label_slival->setObjectName(QString::fromUtf8("label_slival"));

        gridLayout_imginfo->addWidget(label_slival, 1, 0, 1, 1);

        label_slices = new QLabel(gridLayoutWidget);
        label_slices->setObjectName(QString::fromUtf8("label_slices"));

        gridLayout_imginfo->addWidget(label_slices, 0, 0, 1, 1);

        label_acq = new QLabel(gridLayoutWidget);
        label_acq->setObjectName(QString::fromUtf8("label_acq"));

        gridLayout_imginfo->addWidget(label_acq, 0, 1, 1, 1);

        label_acqval = new QLabel(gridLayoutWidget);
        label_acqval->setObjectName(QString::fromUtf8("label_acqval"));

        gridLayout_imginfo->addWidget(label_acqval, 1, 1, 1, 1);

        label_rep = new QLabel(gridLayoutWidget);
        label_rep->setObjectName(QString::fromUtf8("label_rep"));

        gridLayout_imginfo->addWidget(label_rep, 0, 2, 1, 1);

        label_repval = new QLabel(gridLayoutWidget);
        label_repval->setObjectName(QString::fromUtf8("label_repval"));

        gridLayout_imginfo->addWidget(label_repval, 1, 2, 1, 1);

        label_par = new QLabel(gridLayoutWidget);
        label_par->setObjectName(QString::fromUtf8("label_par"));

        gridLayout_imginfo->addWidget(label_par, 0, 3, 1, 1);

        label_parval = new QLabel(gridLayoutWidget);
        label_parval->setObjectName(QString::fromUtf8("label_parval"));

        gridLayout_imginfo->addWidget(label_parval, 1, 3, 1, 1);

        label_echo = new QLabel(gridLayoutWidget);
        label_echo->setObjectName(QString::fromUtf8("label_echo"));

        gridLayout_imginfo->addWidget(label_echo, 0, 4, 1, 1);

        label_echoval = new QLabel(gridLayoutWidget);
        label_echoval->setObjectName(QString::fromUtf8("label_echoval"));

        gridLayout_imginfo->addWidget(label_echoval, 1, 4, 1, 1);

        label_pha = new QLabel(gridLayoutWidget);
        label_pha->setObjectName(QString::fromUtf8("label_pha"));

        gridLayout_imginfo->addWidget(label_pha, 0, 5, 1, 1);

        label_phaval = new QLabel(gridLayoutWidget);
        label_phaval->setObjectName(QString::fromUtf8("label_phaval"));

        gridLayout_imginfo->addWidget(label_phaval, 1, 5, 1, 1);

        label_set = new QLabel(gridLayoutWidget);
        label_set->setObjectName(QString::fromUtf8("label_set"));

        gridLayout_imginfo->addWidget(label_set, 0, 6, 1, 1);

        label_setval = new QLabel(gridLayoutWidget);
        label_setval->setObjectName(QString::fromUtf8("label_setval"));

        gridLayout_imginfo->addWidget(label_setval, 1, 6, 1, 1);

        label_seg = new QLabel(gridLayoutWidget);
        label_seg->setObjectName(QString::fromUtf8("label_seg"));

        gridLayout_imginfo->addWidget(label_seg, 0, 7, 1, 1);

        label_segval = new QLabel(gridLayoutWidget);
        label_segval->setObjectName(QString::fromUtf8("label_segval"));

        gridLayout_imginfo->addWidget(label_segval, 1, 7, 1, 1);


        gridLayout_2->addWidget(frame_4, 0, 1, 2, 1);

        groupBox = new QGroupBox(centralWidget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setGeometry(QRect(10, 120, 631, 141));
        gridLayout = new QGridLayout(groupBox);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(2);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        cb_cropFOV = new QCheckBox(groupBox);
        cb_cropFOV->setObjectName(QString::fromUtf8("cb_cropFOV"));
        QFont font;
        font.setPointSize(10);
        cb_cropFOV->setFont(font);
        cb_cropFOV->setAutoFillBackground(false);
        cb_cropFOV->setCheckable(true);
        cb_cropFOV->setChecked(false);

        verticalLayout_3->addWidget(cb_cropFOV);


        verticalLayout->addLayout(verticalLayout_3);

        radioButton_rawdata = new QRadioButton(groupBox);
        radioButton_rawdata->setObjectName(QString::fromUtf8("radioButton_rawdata"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(radioButton_rawdata->sizePolicy().hasHeightForWidth());
        radioButton_rawdata->setSizePolicy(sizePolicy1);
        radioButton_rawdata->setFont(font);
        radioButton_rawdata->setChecked(true);

        verticalLayout->addWidget(radioButton_rawdata);

        radioButton_ifft = new QRadioButton(groupBox);
        radioButton_ifft->setObjectName(QString::fromUtf8("radioButton_ifft"));
        sizePolicy1.setHeightForWidth(radioButton_ifft->sizePolicy().hasHeightForWidth());
        radioButton_ifft->setSizePolicy(sizePolicy1);
        radioButton_ifft->setFont(font);
        radioButton_ifft->setChecked(false);

        verticalLayout->addWidget(radioButton_ifft);

        radioButton_irgnl2 = new QRadioButton(groupBox);
        radioButton_irgnl2->setObjectName(QString::fromUtf8("radioButton_irgnl2"));
        sizePolicy1.setHeightForWidth(radioButton_irgnl2->sizePolicy().hasHeightForWidth());
        radioButton_irgnl2->setSizePolicy(sizePolicy1);
        radioButton_irgnl2->setFont(font);

        verticalLayout->addWidget(radioButton_irgnl2);

        radioButton_irgntv = new QRadioButton(groupBox);
        radioButton_irgntv->setObjectName(QString::fromUtf8("radioButton_irgntv"));
        sizePolicy1.setHeightForWidth(radioButton_irgntv->sizePolicy().hasHeightForWidth());
        radioButton_irgntv->setSizePolicy(sizePolicy1);
        radioButton_irgntv->setFont(font);

        verticalLayout->addWidget(radioButton_irgntv);

        radioButton_irgntgv = new QRadioButton(groupBox);
        radioButton_irgntgv->setObjectName(QString::fromUtf8("radioButton_irgntgv"));
        sizePolicy1.setHeightForWidth(radioButton_irgntgv->sizePolicy().hasHeightForWidth());
        radioButton_irgntgv->setSizePolicy(sizePolicy1);
        radioButton_irgntgv->setFont(font);

        verticalLayout->addWidget(radioButton_irgntgv);


        gridLayout->addLayout(verticalLayout, 1, 0, 1, 1);

        frame = new QFrame(groupBox);
        frame->setObjectName(QString::fromUtf8("frame"));
        frame->setMinimumSize(QSize(39, 0));
        frame->setFrameShape(QFrame::NoFrame);
        frame->setFrameShadow(QFrame::Raised);

        gridLayout->addWidget(frame, 1, 1, 1, 1);

        gridLayout_imgselect = new QGridLayout();
        gridLayout_imgselect->setSpacing(6);
        gridLayout_imgselect->setObjectName(QString::fromUtf8("gridLayout_imgselect"));
        gridLayout_imgselect->setSizeConstraint(QLayout::SetNoConstraint);
        label_slices_2 = new QLabel(groupBox);
        label_slices_2->setObjectName(QString::fromUtf8("label_slices_2"));
        sizePolicy1.setHeightForWidth(label_slices_2->sizePolicy().hasHeightForWidth());
        label_slices_2->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(label_slices_2, 0, 0, 1, 1);

        label_acq_2 = new QLabel(groupBox);
        label_acq_2->setObjectName(QString::fromUtf8("label_acq_2"));
        sizePolicy1.setHeightForWidth(label_acq_2->sizePolicy().hasHeightForWidth());
        label_acq_2->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(label_acq_2, 0, 1, 1, 1);

        label_rep_2 = new QLabel(groupBox);
        label_rep_2->setObjectName(QString::fromUtf8("label_rep_2"));
        sizePolicy1.setHeightForWidth(label_rep_2->sizePolicy().hasHeightForWidth());
        label_rep_2->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(label_rep_2, 0, 2, 1, 1);

        label_echo_2 = new QLabel(groupBox);
        label_echo_2->setObjectName(QString::fromUtf8("label_echo_2"));
        sizePolicy1.setHeightForWidth(label_echo_2->sizePolicy().hasHeightForWidth());
        label_echo_2->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(label_echo_2, 0, 4, 1, 1);

        label_pha_2 = new QLabel(groupBox);
        label_pha_2->setObjectName(QString::fromUtf8("label_pha_2"));
        sizePolicy1.setHeightForWidth(label_pha_2->sizePolicy().hasHeightForWidth());
        label_pha_2->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(label_pha_2, 0, 5, 1, 1);

        label_set_2 = new QLabel(groupBox);
        label_set_2->setObjectName(QString::fromUtf8("label_set_2"));
        sizePolicy1.setHeightForWidth(label_set_2->sizePolicy().hasHeightForWidth());
        label_set_2->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(label_set_2, 0, 6, 1, 1);

        label_seg_2 = new QLabel(groupBox);
        label_seg_2->setObjectName(QString::fromUtf8("label_seg_2"));
        sizePolicy1.setHeightForWidth(label_seg_2->sizePolicy().hasHeightForWidth());
        label_seg_2->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(label_seg_2, 0, 7, 1, 1);

        combo_sli_select = new QComboBox(groupBox);
        combo_sli_select->setObjectName(QString::fromUtf8("combo_sli_select"));
        sizePolicy1.setHeightForWidth(combo_sli_select->sizePolicy().hasHeightForWidth());
        combo_sli_select->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(combo_sli_select, 1, 0, 1, 1);

        combo_acq_select = new QComboBox(groupBox);
        combo_acq_select->setObjectName(QString::fromUtf8("combo_acq_select"));
        sizePolicy1.setHeightForWidth(combo_acq_select->sizePolicy().hasHeightForWidth());
        combo_acq_select->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(combo_acq_select, 1, 1, 1, 1);

        combo_rep_select = new QComboBox(groupBox);
        combo_rep_select->setObjectName(QString::fromUtf8("combo_rep_select"));
        sizePolicy1.setHeightForWidth(combo_rep_select->sizePolicy().hasHeightForWidth());
        combo_rep_select->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(combo_rep_select, 1, 2, 1, 1);

        combo_par_select = new QComboBox(groupBox);
        combo_par_select->setObjectName(QString::fromUtf8("combo_par_select"));
        sizePolicy1.setHeightForWidth(combo_par_select->sizePolicy().hasHeightForWidth());
        combo_par_select->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(combo_par_select, 1, 3, 1, 1);

        combo_echo_select = new QComboBox(groupBox);
        combo_echo_select->setObjectName(QString::fromUtf8("combo_echo_select"));
        sizePolicy1.setHeightForWidth(combo_echo_select->sizePolicy().hasHeightForWidth());
        combo_echo_select->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(combo_echo_select, 1, 4, 1, 1);

        combo_pha_select = new QComboBox(groupBox);
        combo_pha_select->setObjectName(QString::fromUtf8("combo_pha_select"));
        sizePolicy1.setHeightForWidth(combo_pha_select->sizePolicy().hasHeightForWidth());
        combo_pha_select->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(combo_pha_select, 1, 5, 1, 1);

        combo_seg_select = new QComboBox(groupBox);
        combo_seg_select->setObjectName(QString::fromUtf8("combo_seg_select"));
        sizePolicy1.setHeightForWidth(combo_seg_select->sizePolicy().hasHeightForWidth());
        combo_seg_select->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(combo_seg_select, 1, 7, 1, 1);

        combo_set_select = new QComboBox(groupBox);
        combo_set_select->setObjectName(QString::fromUtf8("combo_set_select"));
        sizePolicy1.setHeightForWidth(combo_set_select->sizePolicy().hasHeightForWidth());
        combo_set_select->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(combo_set_select, 1, 6, 1, 1);

        label_par_2 = new QLabel(groupBox);
        label_par_2->setObjectName(QString::fromUtf8("label_par_2"));
        sizePolicy1.setHeightForWidth(label_par_2->sizePolicy().hasHeightForWidth());
        label_par_2->setSizePolicy(sizePolicy1);

        gridLayout_imgselect->addWidget(label_par_2, 0, 3, 1, 1);


        gridLayout->addLayout(gridLayout_imgselect, 1, 2, 1, 1);

        frame_2 = new QFrame(groupBox);
        frame_2->setObjectName(QString::fromUtf8("frame_2"));
        frame_2->setMinimumSize(QSize(14, 0));
        frame_2->setFrameShape(QFrame::NoFrame);
        frame_2->setFrameShadow(QFrame::Raised);

        gridLayout->addWidget(frame_2, 1, 3, 1, 1);

        layoutWidget = new QWidget(centralWidget);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(10, 260, 621, 231));
        horizontalLayout = new QHBoxLayout(layoutWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        groupBox_3 = new QGroupBox(layoutWidget);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        verticalLayoutWidget_2 = new QWidget(groupBox_3);
        verticalLayoutWidget_2->setObjectName(QString::fromUtf8("verticalLayoutWidget_2"));
        verticalLayoutWidget_2->setGeometry(QRect(0, 20, 301, 211));
        verticalLayout_2 = new QVBoxLayout(verticalLayoutWidget_2);
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setContentsMargins(11, 11, 11, 11);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        text_info = new QTextBrowser(verticalLayoutWidget_2);
        text_info->setObjectName(QString::fromUtf8("text_info"));

        verticalLayout_2->addWidget(text_info);


        horizontalLayout->addWidget(groupBox_3);

        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        verticalSpacer = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout_4->addItem(verticalSpacer);

        groupBox_5 = new QGroupBox(layoutWidget);
        groupBox_5->setObjectName(QString::fromUtf8("groupBox_5"));
        formLayoutWidget = new QWidget(groupBox_5);
        formLayoutWidget->setObjectName(QString::fromUtf8("formLayoutWidget"));
        formLayoutWidget->setGeometry(QRect(70, 0, 234, 111));
        formLayout = new QFormLayout(formLayoutWidget);
        formLayout->setSpacing(6);
        formLayout->setContentsMargins(11, 11, 11, 11);
        formLayout->setObjectName(QString::fromUtf8("formLayout"));
        formLayout->setHorizontalSpacing(2);
        formLayout->setVerticalSpacing(4);
        formLayout->setContentsMargins(6, 0, 0, 0);
        pb_outputpath = new QPushButton(formLayoutWidget);
        pb_outputpath->setObjectName(QString::fromUtf8("pb_outputpath"));

        formLayout->setWidget(4, QFormLayout::FieldRole, pb_outputpath);

        pb_automatic = new QPushButton(formLayoutWidget);
        pb_automatic->setObjectName(QString::fromUtf8("pb_automatic"));
        pb_automatic->setCheckable(false);

        formLayout->setWidget(2, QFormLayout::FieldRole, pb_automatic);

        label = new QLabel(formLayoutWidget);
        label->setObjectName(QString::fromUtf8("label"));

        formLayout->setWidget(2, QFormLayout::LabelRole, label);

        label_2 = new QLabel(formLayoutWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        formLayout->setWidget(4, QFormLayout::LabelRole, label_2);

        cb_savingtype = new QComboBox(formLayoutWidget);
        cb_savingtype->setObjectName(QString::fromUtf8("cb_savingtype"));

        formLayout->setWidget(3, QFormLayout::FieldRole, cb_savingtype);

        label_3 = new QLabel(formLayoutWidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        formLayout->setWidget(3, QFormLayout::LabelRole, label_3);


        verticalLayout_4->addWidget(groupBox_5);

        groupBox_4 = new QGroupBox(layoutWidget);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(groupBox_4->sizePolicy().hasHeightForWidth());
        groupBox_4->setSizePolicy(sizePolicy2);
        formLayoutWidget_3 = new QWidget(groupBox_4);
        formLayoutWidget_3->setObjectName(QString::fromUtf8("formLayoutWidget_3"));
        formLayoutWidget_3->setGeometry(QRect(70, 40, 231, 60));
        formLayout_2 = new QFormLayout(formLayoutWidget_3);
        formLayout_2->setSpacing(6);
        formLayout_2->setContentsMargins(11, 11, 11, 11);
        formLayout_2->setObjectName(QString::fromUtf8("formLayout_2"));
        formLayout_2->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
        formLayout_2->setContentsMargins(0, 0, 0, 0);
        pushButton = new QPushButton(formLayoutWidget_3);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));

        formLayout_2->setWidget(0, QFormLayout::LabelRole, pushButton);

        cb_autofilename = new QCheckBox(formLayoutWidget_3);
        cb_autofilename->setObjectName(QString::fromUtf8("cb_autofilename"));
        cb_autofilename->setChecked(true);

        formLayout_2->setWidget(0, QFormLayout::FieldRole, cb_autofilename);

        progressBar = new QProgressBar(formLayoutWidget_3);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        QSizePolicy sizePolicy3(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(progressBar->sizePolicy().hasHeightForWidth());
        progressBar->setSizePolicy(sizePolicy3);
        progressBar->setLayoutDirection(Qt::LeftToRight);
        progressBar->setAutoFillBackground(false);
        progressBar->setMinimum(0);
        progressBar->setMaximum(100);
        progressBar->setValue(0);
        progressBar->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        progressBar->setTextVisible(false);

        formLayout_2->setWidget(1, QFormLayout::SpanningRole, progressBar);


        verticalLayout_4->addWidget(groupBox_4);


        horizontalLayout->addLayout(verticalLayout_4);

        IMT_MRGui->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(IMT_MRGui);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 687, 23));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuIRGN_Config = new QMenu(menuBar);
        menuIRGN_Config->setObjectName(QString::fromUtf8("menuIRGN_Config"));
        IMT_MRGui->setMenuBar(menuBar);
        mainToolBar = new QToolBar(IMT_MRGui);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        IMT_MRGui->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(IMT_MRGui);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        IMT_MRGui->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuIRGN_Config->menuAction());
        menuFile->addAction(actionRead_DicomRaw);
        menuFile->addAction(actionOpen_meas_dat);
        menuFile->addAction(actionOpen_matlab_bin);
        menuIRGN_Config->addAction(actionParameter);
        menuIRGN_Config->addSeparator();
        menuIRGN_Config->addAction(actionAuto_Archive);
        menuIRGN_Config->addSeparator();
        menuIRGN_Config->addAction(actionSpecial);
        mainToolBar->addAction(actionRead_DicomRaw);
        mainToolBar->addAction(actionOpen_meas_dat);
        mainToolBar->addAction(actionOpen_matlab_bin);

        retranslateUi(IMT_MRGui);

        QMetaObject::connectSlotsByName(IMT_MRGui);
    } // setupUi

    void retranslateUi(QMainWindow *IMT_MRGui)
    {
        IMT_MRGui->setWindowTitle(QApplication::translate("IMT_MRGui", "IMT_MRGui", 0, QApplication::UnicodeUTF8));
        actionRead_DicomRaw->setText(QApplication::translate("IMT_MRGui", "Read DicomRaw", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionRead_DicomRaw->setToolTip(QApplication::translate("IMT_MRGui", "Read DicomRaw [Ctrl+D]", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionRead_DicomRaw->setShortcut(QApplication::translate("IMT_MRGui", "Ctrl+D", 0, QApplication::UnicodeUTF8));
        actionOpen_meas_dat->setText(QApplication::translate("IMT_MRGui", "Open meas.dat", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionOpen_meas_dat->setToolTip(QApplication::translate("IMT_MRGui", "Open meas.dat [Ctrl+M]", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionOpen_meas_dat->setShortcut(QApplication::translate("IMT_MRGui", "Ctrl+M", 0, QApplication::UnicodeUTF8));
        actionParameter->setText(QApplication::translate("IMT_MRGui", "IRGN Parameter", 0, QApplication::UnicodeUTF8));
        actionOpen_matlab_bin->setText(QApplication::translate("IMT_MRGui", "Open matlab.bin", 0, QApplication::UnicodeUTF8));
        actionOpen_matlab_bin->setShortcut(QApplication::translate("IMT_MRGui", "Ctrl+B", 0, QApplication::UnicodeUTF8));
        actionAuto_Archive->setText(QApplication::translate("IMT_MRGui", "Path Settings", 0, QApplication::UnicodeUTF8));
        actionSpecial->setText(QApplication::translate("IMT_MRGui", "special", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("IMT_MRGui", "RawData Information", 0, QApplication::UnicodeUTF8));
        label_rows->setText(QApplication::translate("IMT_MRGui", "Rows: ", 0, QApplication::UnicodeUTF8));
        label_columns->setText(QApplication::translate("IMT_MRGui", "Columns:", 0, QApplication::UnicodeUTF8));
        label_coils->setText(QApplication::translate("IMT_MRGui", "Coils:", 0, QApplication::UnicodeUTF8));
        label_rowsval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
        label_columnsval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
        label_coilsval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
        label_slival->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_slices->setToolTip(QApplication::translate("IMT_MRGui", "Number of Slices", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_slices->setText(QApplication::translate("IMT_MRGui", "Sli", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_acq->setToolTip(QApplication::translate("IMT_MRGui", "Number of Acquisitions", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_acq->setText(QApplication::translate("IMT_MRGui", "Acq", 0, QApplication::UnicodeUTF8));
        label_acqval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_rep->setToolTip(QApplication::translate("IMT_MRGui", "Number of Repititions", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_rep->setText(QApplication::translate("IMT_MRGui", "Rep", 0, QApplication::UnicodeUTF8));
        label_repval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_par->setToolTip(QApplication::translate("IMT_MRGui", "Number of Partitions", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_par->setText(QApplication::translate("IMT_MRGui", "Par", 0, QApplication::UnicodeUTF8));
        label_parval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_echo->setToolTip(QApplication::translate("IMT_MRGui", "Number of Echos", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_echo->setText(QApplication::translate("IMT_MRGui", "Echo", 0, QApplication::UnicodeUTF8));
        label_echoval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_pha->setToolTip(QApplication::translate("IMT_MRGui", "Number of Phases", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_pha->setText(QApplication::translate("IMT_MRGui", "Pha", 0, QApplication::UnicodeUTF8));
        label_phaval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_set->setToolTip(QApplication::translate("IMT_MRGui", "Number of Sets", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_set->setText(QApplication::translate("IMT_MRGui", "Set", 0, QApplication::UnicodeUTF8));
        label_setval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_seg->setToolTip(QApplication::translate("IMT_MRGui", "Number of Segments", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_seg->setText(QApplication::translate("IMT_MRGui", "Seg", 0, QApplication::UnicodeUTF8));
        label_segval->setText(QApplication::translate("IMT_MRGui", "---", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("IMT_MRGui", "Configuration", 0, QApplication::UnicodeUTF8));
        cb_cropFOV->setText(QApplication::translate("IMT_MRGui", "Crop FOV", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        radioButton_rawdata->setToolTip(QApplication::translate("IMT_MRGui", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Ubuntu'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">No Calculation</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Rawdata will be saved</p></body></html>", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        radioButton_rawdata->setText(QApplication::translate("IMT_MRGui", "Rawdata", 0, QApplication::UnicodeUTF8));
        radioButton_ifft->setText(QApplication::translate("IMT_MRGui", "Inverse FFT", 0, QApplication::UnicodeUTF8));
        radioButton_irgnl2->setText(QApplication::translate("IMT_MRGui", "IRGN L2", 0, QApplication::UnicodeUTF8));
        radioButton_irgntv->setText(QApplication::translate("IMT_MRGui", "IRGN TV", 0, QApplication::UnicodeUTF8));
        radioButton_irgntgv->setText(QApplication::translate("IMT_MRGui", "IRGN TGV", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_slices_2->setToolTip(QApplication::translate("IMT_MRGui", "Select Slice", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_slices_2->setText(QApplication::translate("IMT_MRGui", " Sli", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_acq_2->setToolTip(QApplication::translate("IMT_MRGui", "Select Acquisition", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_acq_2->setText(QApplication::translate("IMT_MRGui", "Acq", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_rep_2->setToolTip(QApplication::translate("IMT_MRGui", "Select Repitition", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_rep_2->setText(QApplication::translate("IMT_MRGui", "Rep", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_echo_2->setToolTip(QApplication::translate("IMT_MRGui", "Select Echo", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_echo_2->setText(QApplication::translate("IMT_MRGui", "Echo", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_pha_2->setToolTip(QApplication::translate("IMT_MRGui", "Select Phase", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_pha_2->setText(QApplication::translate("IMT_MRGui", "Pha", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_set_2->setToolTip(QApplication::translate("IMT_MRGui", "Select Set", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_set_2->setText(QApplication::translate("IMT_MRGui", "Set", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_seg_2->setToolTip(QApplication::translate("IMT_MRGui", "Select Segment", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_seg_2->setText(QApplication::translate("IMT_MRGui", "Seg", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        label_par_2->setToolTip(QApplication::translate("IMT_MRGui", "Select Partition", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_par_2->setText(QApplication::translate("IMT_MRGui", "Par", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("IMT_MRGui", "Info", 0, QApplication::UnicodeUTF8));
        groupBox_5->setTitle(QString());
        pb_outputpath->setText(QApplication::translate("IMT_MRGui", "Output Path", 0, QApplication::UnicodeUTF8));
        pb_automatic->setText(QApplication::translate("IMT_MRGui", "Automatic", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("IMT_MRGui", "Operation mode:", 0, QApplication::UnicodeUTF8));
        label_2->setText(QString());
        cb_savingtype->clear();
        cb_savingtype->insertItems(0, QStringList()
         << QApplication::translate("IMT_MRGui", "real", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("IMT_MRGui", "imag", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("IMT_MRGui", "abs", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("IMT_MRGui", "phase", 0, QApplication::UnicodeUTF8)
        );
#ifndef QT_NO_TOOLTIP
        cb_savingtype->setToolTip(QApplication::translate("IMT_MRGui", "Saving Options:\n"
" abs        Absolute Value\n"
" phase   Phase Value\n"
" imag     Imaginary Value\n"
" real       Real Value", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label_3->setText(QApplication::translate("IMT_MRGui", "Saving Type:", 0, QApplication::UnicodeUTF8));
        groupBox_4->setTitle(QString());
        pushButton->setText(QApplication::translate("IMT_MRGui", "Start", 0, QApplication::UnicodeUTF8));
        cb_autofilename->setText(QApplication::translate("IMT_MRGui", "Auto Filename", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("IMT_MRGui", "File", 0, QApplication::UnicodeUTF8));
        menuIRGN_Config->setTitle(QApplication::translate("IMT_MRGui", "Settings", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class IMT_MRGui: public Ui_IMT_MRGui {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IMT_MRGUI_H
