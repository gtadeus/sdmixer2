#include "sdmixer.h"
#include "ui_sdmixer.h"

#include <QFileDialog>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>

#include <QInputDialog>
#include <QMimeData>
#include <QDragEnterEvent>
#include <QDebug>

#include <fcntl.h>
#include <unistd.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>

sdmixer::sdmixer(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::sdmixer)
{
    ui->setupUi(this);

    setAcceptDrops(true);
    listWidget=ui->listWidget;

}

sdmixer::~sdmixer()
{
    delete ui;
}
void sdmixer::dragEnterEvent(QDragEnterEvent *e)
{
    if (e->mimeData()->hasUrls()) {
        e->acceptProposedAction();
    }
}

void sdmixer::dropEvent(QDropEvent *e)
{
    foreach (const QUrl &url, e->mimeData()->urls()) {
        const QString &fileName = url.toLocalFile();
        qDebug() << "Dropped file:" << fileName;
        insertItem(fileName);
    }
}
void sdmixer::insertItem(QString filename)
{
    QListWidgetItem *item = new QListWidgetItem(QIcon(filename), filename, listWidget);

    /*QString itemText = QInputDialog::getText(this, tr("Insert Item"),
        tr("Input text for the new item:"));

    if (itemText.isNull())
        return;


    QListWidgetItem *newItem = new QListWidgetItem;
    newItem->setText(itemText);

    int row = listWidget->row(listWidget->currentItem());

    listWidget->insertItem(row, newItem);*/

}


void sdmixer::on_actionOpen_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QString(),
            tr("Localisations File (*.txt);;PairFinder Files (*.out)"));

    if (!fileName.isEmpty())
    {
        int lines = read_file(fileName.toLatin1());
        if ( lines != 0)
        {
            QMessageBox::information(this, tr("Information"), QString::number(lines));
        }else
        {
            QMessageBox::critical(this, tr("Error"), tr("Could not open file"));
        }
    }

}

int sdmixer::read_file(char const *fname)
{
    int lines = 0;
    int cols = 0;
    double dd;

    static const int BUFFER_SIZE = 16*1024;
    int fd = open(fname, O_RDONLY);
    if(fd == -1)
        return 0;
    char buf[BUFFER_SIZE + 1];
    while(size_t bytes_read = read(fd, buf, BUFFER_SIZE))
    {
        if(bytes_read == (size_t)-1)
            return 0;//handle_error("read failed");
        if (!bytes_read)
            break;

        for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p)
           ++lines;
    }

    std::ifstream ifs(fname);
    std::string firstLine, secondLine;
    getline (ifs, firstLine);
    getline (ifs, secondLine);
    ifs.close();

    std::stringstream countColsStream(secondLine);

    while (countColsStream >> dd)
    {
        ++cols;
    }


    FILE *infile = fopen(fname, "r");
    char buffer[4096];

    gsl_matrix *K = gsl_matrix_alloc(lines,cols);
   // double* a = new double[lines*cols];

    int curr_line=0;

    fseek ( infile , int(firstLine.length()), SEEK_SET );

    while (fgets(buffer, sizeof(buffer), infile))
    {
            double d;
            std::stringstream lineStream(buffer);

            int curr_col = 0;
            while (lineStream >> d)
            {
                //ui->textEdit->append(QString::number(d) + "  " +  QString::number(curr_line) +  "  " + QString::number(curr_col) + "\n" );
                gsl_matrix_set(K, curr_line, curr_col, d);
                //a[curr_line*cols+curr_col]=d;
                ++curr_col;
            }

            ++curr_line;

    }


    QString qs(firstLine.c_str());

    ui->textEdit->append(QString::number(lines));
    ui->textEdit->append(QString::number(cols));
    ui->textEdit->append(qs);
    return lines;
}

void sdmixer::on_actionQuit_triggered()
{
    qApp->quit();
}



void sdmixer::on_addFileButton_clicked()
{

}

