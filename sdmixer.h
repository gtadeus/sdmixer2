#ifndef SDMIXER_H
#define SDMIXER_H

#define DEFAULT_SETTINGS "default_settings.txt"

#include <QMainWindow>
#include <QListWidget>
#include <QDebug>
#include <string>
#include <vector>
#include <QString>
#include <QTextEdit>
#include <QTextStream>
#include <ctime>

class Settings;

namespace Ui {
class sdmixer;
}

class sdmixer : public QMainWindow
{
    Q_OBJECT

public:
    struct fishing{
        bool run=false;
        int increment=0;
        int range=0;
        int subset=0;
    };

    explicit sdmixer(QWidget *parent = 0);
    ~sdmixer();
    int read_file(const char *file);
    void insertItem(QString filename, QListWidget *list);
    QString timestamp();
    QString error_msg(QString msg);
    void writeToConsole(QString q);
    bool getSettingsFromUI();
    void setSettingsToUI(Settings s);
    bool prepareForRun();


    bool getRunPairfinder();
    bool getRunFilter();
    bool getRunReconstructor();
    bool getForce2D();
    int getPixelSize();
    double getOffset(int dim);
    double getEpsilon(int dim);
    fishing getFishing();
    std::vector<QString> getInputFiles();
    std::vector<QString> getFilterFiles();
    double getMaxIntShort();
    double getMaxIntLong();
    double getPrecision();
    double getReconstructor_xyBinning();
    double getReconstructor_zBinning();

protected:
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);
private slots:
    void on_actionQuit_triggered();
    void on_addFileButton_clicked();
    void on_startDemixing_clicked();
    void on_actionAbout_sdmixer_triggered();

    void on_removeFilterButton_clicked();

    void on_removeFilesButton_clicked();

    void on_actionSave_Settings_as_Default_triggered();

    void on_actionLoad_Settings_triggered();

    void on_actionSave_Preferences_triggered();

private:

    std::vector<QString> InputFiles;
    bool runPairFinder=false;
    bool runFilter=false;
    bool runReconstructor=false;
    bool force2D=false;
    int pixelSizeNM=0;

    static const int max_dims = 3;
    double offset[max_dims]={0};
    double epsilon[max_dims]={0};

    fishing fishingSettings;
    int maxIntensityLong=0;
    int maxIntensityShort=0;
    double precision=0;
    double xyBinning=0;
    double zBinning=0;


    std::vector<QString> FilterFiles;

    Ui::sdmixer *ui;
    QListWidget *listWidgetInputFiles;
    QListWidget *listWidgetFilters;
    QTextEdit *console;
};

#endif // SDMIXER_H
