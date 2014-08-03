#include "sdmixer.h"
#include "pairfinder.h"
#include "reconstructor.h"
#include "filter.h"
#include "ui_sdmixer.h"
#include "settings.h"

#include <QFileDialog>
#include <QFile>
#include <QMessageBox>

#include <QInputDialog>
#include <QMimeData>
#include <QDragEnterEvent>

#include <QThread>



#include <unistd.h>





sdmixer::sdmixer(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::sdmixer)
{

    ui->setupUi(this);
    setAcceptDrops(true);

    listWidgetInputFiles = ui->listWidget_InputFiles;
    listWidgetFilters = ui->listWidget_FilterFiles;
    console = ui->textConsole;
    console->acceptRichText();

    console->append(timestamp() + " : Started sdmixer");
    //console->append(error_msg("some error msg"));

    QFile file(DEFAULT_SETTINGS);
    if (file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qDebug()<< "loading default settings...";
        file.close();
        Settings s(DEFAULT_SETTINGS);
        setSettingsToUI(s);

    }
    else
    {
        qDebug() << "no default settings available";
    }

}
//small helper functions
QString sdmixer::timestamp()
{
    QString retval;
    QTextStream out(&retval);

    std::time_t t = std::time(NULL);
    char c_string[100];
    if (std::strftime(c_string, sizeof(c_string), "%Y-%m-%d %X", std::localtime(&t)))
    {
        out << c_string;
    }
    return retval;
}
QString sdmixer::error_msg(QString msg)
{
    QString retval;
    QTextStream out(&retval);
    out << "<font color=\"#FF0000\">" << timestamp() << " : ERROR! " << msg << "</font><br>";
    return retval;
}

void sdmixer::writeToConsole(QString q)
{
    console->append(q);
}

sdmixer::~sdmixer()
{
    delete ui;
}

//Get everything from ui

bool sdmixer::getSettingsFromUI()
{
    InputFiles.clear();

    for(int i = 0; i < ui->listWidget_InputFiles->count(); ++i)
    {
        QListWidgetItem* item = ui->listWidget_InputFiles->item(i);
        InputFiles.push_back(item->text());
    }

    runPairFinder = ui->checkBox_runPairFinder->isChecked();
    runFilter = ui->checkBox_runFilter->isChecked();
    runReconstructor = ui->checkBox_runReconstructor->isChecked();

    force2D = ui->checkBox_force2D->isChecked();

    pixelSizeNM = QString(ui->lineEdit_pixelSizeNM->text()).toInt();

    QString temp;
    temp = ui->lineEdit_xOffset->text();
    if( !temp.isEmpty() )
    {
        offset[0]=temp.toDouble();
        temp = ui->lineEdit_yOffset->text();
        if( ! temp.isEmpty() )
        {
            offset[1]=temp.toDouble();
            temp = ui->lineEdit_zOffset->text();

            if( ! temp.isEmpty() )
            {
                offset[2]=temp.toDouble();
            }
        }
    }

    temp = ui->lineEdit_xEpsilon->text();
    if( !temp.isEmpty() )
    {
        epsilon[0]=temp.toDouble();
        temp = ui->lineEdit_yEpsilon->text();
        if( ! temp.isEmpty() )
        {
            epsilon[1]=temp.toDouble();
            temp = ui->lineEdit_zEpsilon->text();

            if( ! temp.isEmpty() )
            {
                epsilon[2]=temp.toDouble();
            }
        }
    }

    fishingSettings.run = ui->checkBox_runFishing->isChecked() ;
    if (fishingSettings.run)
    {
        fishingSettings.increment = QString(ui->lineEdit_FishingIncrement->text()).toInt();
        fishingSettings.range = QString(ui->lineEdit_FishingRange->text()).toInt();
        fishingSettings.subset = QString(ui->lineEdit_FishingSubset->text()).toInt();
    }


    maxIntensityLong =  QString(ui->lineEdit_maxIntLong->text()).toDouble();
    maxIntensityShort = QString(ui->lineEdit_maxIntShort->text()).toDouble();
    precision = QString(ui->lineEdit_precision->text()).toDouble();

    xyBinning = QString(ui->lineEdit_xyBinning->text()).toDouble();
    zBinning = QString(ui->lineEdit_zBinning->text()).toDouble();

    FilterFiles.clear();
    for(int i = 0; i < ui->listWidget_FilterFiles->count(); ++i)
    {
        QListWidgetItem* item = ui->listWidget_FilterFiles->item(i);
        FilterFiles.push_back(item->text());
    }

    return true;

}
void sdmixer::setSettingsToUI(Settings s){
    ui->lineEdit_pixelSizeNM->setText(QString::number(s.getPixelSizeNM()));
    ui->checkBox_runPairFinder->setChecked(s.getRunPairFinder());
    ui->checkBox_runFilter->setChecked(s.getRunFilter());
    ui->checkBox_runPairFinder->setChecked(s.getRunReconstructor());
    ui->checkBox_force2D->setChecked(s.getForce2D());


    ui->lineEdit_xOffset->setText(QString::number(s.getOffset(0)));
    ui->lineEdit_yOffset->setText(QString::number(s.getOffset(1)));
    ui->lineEdit_zOffset->setText(QString::number(s.getOffset(2)));

    ui->lineEdit_xEpsilon->setText(QString::number(s.getEpsilon(0)));
    ui->lineEdit_yEpsilon->setText(QString::number(s.getEpsilon(1)));
    ui->lineEdit_zEpsilon->setText(QString::number(s.getEpsilon(2)));

    fishing f =s.getFishing();

    ui->checkBox_runFishing->setChecked(f.run);
    ui->lineEdit_FishingIncrement->setText(QString::number(f.increment));
    ui->lineEdit_FishingRange->setText(QString::number(f.range));
    ui->lineEdit_FishingSubset->setText(QString::number(f.subset));

    ui->lineEdit_maxIntLong->setText(QString::number(s.getMaxIntLong()));
    ui->lineEdit_maxIntShort->setText(QString::number(s.getMaxIntShort()));
    ui->lineEdit_precision->setText(QString::number(s.getPrecision()));

    ui->lineEdit_xyBinning->setText(QString::number(s.getXYbinning()));
    ui->lineEdit_zBinning->setText(QString::number(s.getZbinning()));
}

void sdmixer::setStartDemixingButtonEnabled(bool val)
{
    ui->startDemixing->setEnabled(val);
}


std::vector<QString> sdmixer::getInputFiles(){return this->InputFiles;}
bool sdmixer::getRunPairfinder(){return this->runPairFinder;}
bool sdmixer::getRunFilter(){return this->runFilter;}
bool sdmixer::getRunReconstructor(){return this->runReconstructor;}
bool sdmixer::getForce2D(){return this->force2D;}
int sdmixer::getPixelSize(){return this->pixelSizeNM;}
double sdmixer::getOffset(int dim){
    qDebug() << "offset " << dim << ": " <<offset[dim];
    return offset[dim];
}
double sdmixer::getEpsilon(int dim){
    qDebug()<<" epsilon " << dim << ": " << epsilon[dim];
    return epsilon[dim];}
sdmixer::fishing sdmixer::getFishing(){return this->fishingSettings;}
std::vector<QString> sdmixer::getFilterFiles(){return this->FilterFiles;}
double sdmixer::getMaxIntShort(){return this->maxIntensityShort;}
double sdmixer::getMaxIntLong(){return this->maxIntensityLong;}
double sdmixer::getPrecision(){return this->precision;}
double sdmixer::getReconstructor_xyBinning(){return this->xyBinning;}
double sdmixer::getReconstructor_zBinning(){return this->zBinning;}




// Drag & Drop implementation
// Input Files can be draged and dropped onto MainWindow
void sdmixer::dragEnterEvent(QDragEnterEvent *e){
    if (e->mimeData()->hasUrls()) {
        e->acceptProposedAction();
    }
}

void sdmixer::dropEvent(QDropEvent *e){
    foreach (const QUrl &url, e->mimeData()->urls()) {
        const QString &fileName = url.toLocalFile();
        qDebug() << "Dropped file:" << fileName;
        insertItem(fileName, listWidgetInputFiles);

    }
}
void sdmixer::insertItem(QString filename, QListWidget *list){
    QListWidgetItem *item = new QListWidgetItem(QIcon(filename), filename, list);
}

void sdmixer::on_actionQuit_triggered(){
    qApp->quit();
}

void sdmixer::on_addFileButton_clicked(){
    QStringList fileName = QFileDialog::getOpenFileNames(this, tr("Open File"), QString(), tr("Localisations File (*.txt);;PairFinder Files (*.out);;All Files (*.*)"));

    if (!fileName.isEmpty())
    {
        for(auto i : fileName)
        {
            insertItem(i, listWidgetInputFiles);
        }
    }
}


void sdmixer::on_startDemixing_clicked()
{
    getSettingsFromUI();

    /*if(InputFiles.empty())
    {
        QString err = error_msg("no input files selected");
        writeToConsole(err);
        return;
    }*/
   // for (auto &i : InputFiles)
    {
        QThread* thread = new QThread;
        PairFinder* p = new PairFinder(this, "random"/*InputFiles[0]*/);

        p->moveToThread(thread);
        connect(thread, SIGNAL(started()), p, SLOT(doWork()));
        connect(p, SIGNAL(finished()), thread, SLOT(quit()));
        connect(p, SIGNAL(finished()), p, SLOT(deleteLater()));
        connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));
        thread->start();

        setStartDemixingButtonEnabled(false);


        qDebug() << p->getDimensions();

    }
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

void sdmixer::on_removeFilterButton_clicked()
{
    qDeleteAll(listWidgetFilters->selectedItems());
}

void sdmixer::on_removeFilesButton_clicked()
{
    qDeleteAll(listWidgetInputFiles->selectedItems());
}

void sdmixer::on_actionSave_Settings_as_Default_triggered()
{
    getSettingsFromUI();
    Settings s(this);
    s.initXML();
    s.writeSettingsToFile(DEFAULT_SETTINGS);
}

void sdmixer::on_actionLoad_Settings_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QString(), tr("Settings Files (*.txt);;All Files (*.*)"));

    if (!fileName.isEmpty())
    {
        Settings s(fileName);
        setSettingsToUI(s);
    }
}

void sdmixer::on_actionSave_Preferences_triggered()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Save Settings", QString(), "All Files (*.*)");
    if(!fileName.isEmpty())
    {
        getSettingsFromUI();
        Settings s(this);
        s.initXML();
        s.writeSettingsToFile(fileName);
    }
}

void sdmixer::on_addFilterButton_clicked()
{
    QStringList fileName = QFileDialog::getOpenFileNames(this, tr("Open File"), QString(), tr("Image Files (*.png);;All Files (*.*)"));

    if (!fileName.isEmpty())
    {
        for(auto i : fileName)
        {
            insertItem(i, listWidgetFilters);
        }
    }

}
