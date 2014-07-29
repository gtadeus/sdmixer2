#ifndef SDMIXER_H
#define SDMIXER_H

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


    bool getRunPairfinder();
    bool getRunFilter();
    bool getRunReconstructor();
    bool getForce2D();
    int getPixelSize();
    std::vector<double> getOffset();
    std::vector<double> getEpsilon();
    fishing getFishing();
    std::vector<QString> getInputFiles();
    std::vector<QString> getFilterFiles();
    double getMaxIntShort();
    double getMaxIntLong();
    double getPrecision();
    double getReconstructor_xyBinning();
    double getReconstructor_zBinning();

    void setInputFiles(std::vector<QString> v);
    void setRunPairfinder(bool val);
    void setRunFilter(bool val);
    void setRunReconstructor(bool val);
    void setForce2D(bool val);
    void setPixelSize(int val);
    void setOffset(std::vector<double> v);
    void setEpsilon(std::vector<double> v);
    void setFishing(fishing f);
    void setFilterFiles(std::vector<QString> v);
    void setMaxIntShort(int val);
    void setMaxIntLong(int val);
    void setPrecision(double val);
    void setReconstructor_xyBinning(double val);
    void setReconstructor_zBinning(double val);

protected:
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);
private slots:

    void on_actionOpen_triggered();
    void on_actionQuit_triggered();
    void on_addFileButton_clicked();
    void on_startDemixing_clicked();
    void on_actionAbout_sdmixer_triggered();

    void on_removeFilterButton_clicked();

    void on_removeFilesButton_clicked();

private:
    std::vector<QString> InputFiles;
    bool runPairFinder=false;
    bool runFilter=false;
    bool runReconstructor=false;
    bool force2D=false;
    int pixelSizeNM=0;
    std::vector<double> offset;
    std::vector<double> epsilon;
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
