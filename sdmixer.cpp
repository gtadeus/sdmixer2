#include "sdmixer.h"
#include "reconstructor.h"
#include "filter.h"
#include "xmlparser.h"
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

    listWidgetInputFiles = ui->listInputFiles;
    listWidgetFilters = ui->listFilters;

}

sdmixer::~sdmixer()
{
    delete ui;
}

//Get everything from ui

bool sdmixer::getRunPairfinder()
{
    return ui->runPairFinder_CheckBox->isChecked();
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
    QListWidgetItem *item = new QListWidgetItem(QIcon(filename), filename, listWidgetInputFiles);

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

    ui->textConsole->append(QString::number(lines));
    ui->textConsole->append(QString::number(cols));
    ui->textConsole->append(qs);
    return lines;
}

void sdmixer::on_actionQuit_triggered()
{
    qApp->quit();
}



void sdmixer::on_addFileButton_clicked()
{
    XMLParser x;
}


void sdmixer::on_startDemixing_clicked()
{
   /* Reconstructor r;
    r.run();*/
    Filter f;
    f.run();

}


void sdmixer::on_actionAbout_sdmixer_triggered()
{
    QMessageBox msgBox(this);
    msgBox.setIcon(QMessageBox::Information);
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setDefaultButton(QMessageBox::Ok);
    msgBox.setText("<p>sdmixer - Analysis of 2D/3D multicolor SD-dSTORM data</p><p>written by Georgi Tadeus at FMP Berlin,<br>Department for Molecular Pharmacology and Cell Biology</p>Feedback is highly appreciated! <a href='mailto:georgi.tadeus@gmail.com?Subject=sdmixer'>georgi.tadeus@gmail.com</a><p>Many thanks to J. Schmoranzer and A. Lampe!</p><p>If you find this tool useful, please cite:</p><p>Lampe, A., Haucke, V., Sigrist, S. J., Heilemann, M. and Schmoranzer, J. (2012), Multi-colour direct STORM with red emitting carbocyanines. Biology of the Cell, 104: 229–237</p>");
    int ret = msgBox.exec();
}
