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
protected:
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);
private slots:

    void on_actionOpen_triggered();

    void on_actionQuit_triggered();

    void on_addFileButton_clicked();

private:
    Ui::sdmixer *ui;
    QListWidget *listWidget;
};

#endif // SDMIXER_H
