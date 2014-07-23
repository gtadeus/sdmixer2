#ifndef SDMIXER_H
#define SDMIXER_H

#include <QMainWindow>
#include <QListWidget>

namespace Ui {
class sdmixer;
}

class sdmixer : public QMainWindow
{
    Q_OBJECT

public:
    explicit sdmixer(QWidget *parent = 0);
    ~sdmixer();
    int read_file(const char *file);
    void insertItem(QString filename);

    bool getRunPairfinder();
    bool getRunFilter();
    bool getRunReconstructor();
    int getPixelSize();

    std::vector<double> getOffset();
    std::vector<double> getEpsilon();


protected:
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);
private slots:

    void on_actionOpen_triggered();

    void on_actionQuit_triggered();

    void on_addFileButton_clicked();

    void on_startDemixing_clicked();


    void on_actionAbout_sdmixer_triggered();

private:
    Ui::sdmixer *ui;
    QListWidget *listWidgetInputFiles;
    QListWidget *listWidgetFilters;
};

#endif // SDMIXER_H
