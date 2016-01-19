#ifndef IMT_MRGUI_H
#define IMT_MRGUI_H

#include <QMainWindow>
#include <QtConcurrentRun>
#include <QFuture>
#include <QFutureWatcher>
#include <QMessageBox>
#include "widgets/pathsetting.h"
#include "widgets/specialsetting.h"
#include "widgets/irgn_para.h"
#include "widgets/cycleread.h"
#include <complex>

// -- std - C/C++ - includes
#include <iostream>
#include <vector>
#include <time.h>

#include "agile/io/readSiemensVD11.hpp"
#include "agile/calc/postprocess.hpp"
#include "agile/calc/fft.hpp"

#include "agile/calc/genkspacefov.hpp"

#include "agile/calc/l2solve.hpp"
#include "agile/calc/tvsolve.hpp"
#include "agile/calc/tgvsolve.hpp"

//typedef for calculation - float or double cuda calculation
typedef std::complex<double> TType;

typedef std::complex<float> TType_data;//typedef for data - float or double data format
typedef float TType_image;  //don't change image type - no need for double-type image quality
typedef agile::to_real_type<TType_data>::type TType_real_data;
typedef agile::to_real_type<TType>::type TType_real;

struct INDEX_select{
  unsigned short int acq;
  unsigned short int sli;
  unsigned short int par;
  unsigned short int echo;
  unsigned short int pha;
  unsigned short int rep;
  unsigned short int set;
  unsigned short int seg;
};


namespace Ui {
    class IMT_MRGui;
}

class IMT_MRGui : public QMainWindow
{
    Q_OBJECT

public:
    explicit IMT_MRGui(QWidget *parent = 0);
    ~IMT_MRGui();

protected:
    void closeEvent(QCloseEvent *event);

private slots:
  void on_actionRead_DicomRaw_triggered();

  void on_actionOpen_meas_dat_triggered();

  void on_pushButton_clicked();

  void on_actionParameter_triggered();

  void set_Info(QString infotext, bool progress);
  void set_Warning(QString warningtext);
  void set_SaveFile(QString title, QString extension);

  void on_pb_outputpath_clicked();

  void startauto(QStringList filepath, QStringList filename);

  void on_pb_automatic_clicked();

  void on_actionOpen_matlab_bin_triggered();

  void on_actionAuto_Archive_triggered();

  void on_actionSpecial_triggered();

signals:
    void sendInfo(QString infotext, bool progress);
    void sendWarning(QString warningtext);
    void sendSaveFile(QString title, QString extension);


private:
    agile::GPUEnvironment GPU0_;

    Ui::IMT_MRGui *ui;
    void setStatusBar();
    void setTitleText(QString text);

    bool maybeSave();
    void gui_readSiemens(QString path, bool dicom);
    void gui_readMatlab(QString path);
    int startcalculation(unsigned num_rows,unsigned num_columns,unsigned num_coils,
                                    const std::vector<TType_data>& matrix_read, std::vector<TType_image>& image_data);
    int writeoutfile(unsigned num_rows, unsigned num_columns,unsigned num_coils, const std::vector<TType_real_data>& image_data,
                     const std::vector<TType_data>& matrix_read, QString name_index);
    void readcalcwrite_thread();
    void generate_idx(agile::ReadSiemensVD11* read_siemens, std::vector<INDEX_select>& index_vector);
    void irgncalc(agile::IRGN<TType,TType_real>* IRGNobj, std::vector<TType_image>& image_data);
    void copytoarchive();
    void get_dcmrawdata(QStringList filenames);
    void calc_postprocess(const agile::GPUMatrix<TType>* in_mat, agile::GPUMatrix<TType_real>& out_mat);
    void simu_undersampling(agile::GPUMatrix<TType>* Coil, unsigned num_coils);
    void createstartcalcthread();

    PathSetting* _pathsetting;
    SpecialSetting* _specialsetting;
    IRGN_para *open_irgnpara_window; // Be sure to destroy you window somewhere

    QLabel *_statusBarLabelleft;
    QLabel *_statusBarLabelright;

//    QString _dcmrawpath;
//    QString _measdatpath;
//    QString _outputpath;
//    QString _matlabpath;
//    QString _archivepath;
//    QString _dcmrawpath_start;
//    QString _measdatpath_start;
//    QString _outputpath_start;
//    QString _matlabpath_start;
//    QString _archivepath_start;

    unsigned _write_file_index;

    QString _outfilename;
    QString _infotext;
    QString _strselection;

    QString _str_outfilename;
    QStringList _filenames_to_archive;


    agile::ReadSiemensVD11* read_siemens;
    agile::IRGN<TType,TType_real>* IRGNobj;
    agile::IRGN_Params irgnpara;

    agile::FFT<TType>* fftobj_;
    agile::KSpaceFOV<TType>* kspacefovobj_;
    //std::vector<TType_real> zero_vec_cpu_nxny_;
    std::vector<TType_real> zero_vec_cpu_;

    QFuture<void> *future;
    QFutureWatcher<void> *watcher;

    bool _dicom;

    unsigned _num_rows_calc;
    unsigned _num_columns_calc;

    unsigned _num_rows_data;
    unsigned _num_columns_data;
    unsigned _num_coils_data;

    CycleRead* _cycleread_thread;

    unsigned _act_index;
    bool _automatic_on;

    agile::DICOM _in_dicomfile;

    agile::PostProcess<TType, TType_real>* _postprocess;

    unsigned _factor_x;
    unsigned _factor_y;
    unsigned _ref_lines;

    bool _readmatlab;
    std::vector<TType_data> _matlab_read;

    unsigned _matlab_num_rows;
    unsigned _matlab_num_columns;
    unsigned _matlab_num_coils;

    time_t _seconds;

};

#endif // IMT_MRGUI_H
