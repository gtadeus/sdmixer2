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
    pf_output = new std::vector<Localization>();

    ui->setupUi(this);
    setAcceptDrops(true);
    ui->pushButton_CancelRun->setVisible(false);

    listWidgetInputFiles = ui->listWidget_InputFiles;

    // get Info for selected Input File
    connect(listWidgetInputFiles, SIGNAL(itemClicked(QListWidgetItem*)),
            this, SLOT(InputFileClicked(QListWidgetItem*)));
    // Camera Orientation is changed, ->Fill Top/Bottom Left/Right Combobox
    connect(ui->comboBox_CameraOrientation, SIGNAL(highlighted(QString)),
            this, SLOT(CameraOrientationChanged(QString)));

    connect(ui->lineEdit_maxIntLong, SIGNAL(textChanged(QString)),
            this, SLOT(filterInfo(QString)));
    connect(ui->lineEdit_maxIntShort, SIGNAL(textChanged(QString)),
            this, SLOT(filterInfo(QString)));
    connect(ui->lineEdit_precision, SIGNAL(textChanged(QString)),
            this, SLOT(filterInfo(QString)));
    connect(ui->comboBox_FilterOrientation, SIGNAL(currentIndexChanged(QString)),
            this, SLOT(filterInfo(QString)));

    connect(ui->listWidget_FilterFiles, SIGNAL(currentRowChanged(int)),
            this, SLOT(loadKernel(int)));
    connect(ui->comboBox_ConvolutionChannel, SIGNAL(highlighted(QString)),
            this, SLOT(KernelChannelHighlighted(QString)));
    connect(ui->checkBox_sameKernel, SIGNAL(stateChanged(int)),
            this, SLOT(SameKernelCheckBoxChanged(int)));

    UpdateShortPositionComboBox(ui->comboBox_CameraOrientation->itemText(0));

    listWidgetFilters = ui->listWidget_FilterFiles;
    console = ui->textConsole;
    console->acceptRichText();
    log(this, "Started sdmixer");


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
void sdmixer::KernelChannelHighlighted(QString str)
{
    qDebug() << "KernelChannelHighlighted " << str;

    std::vector<gaussian_kernel>::iterator it;

    bool foundIt=false;

    for( it = vec_kernel.begin(); it != vec_kernel.end(); ++it)
    {
        gaussian_kernel gk = *it;
        if(!gk.filterName.compare(str))
        {
            foundIt=true;
            ui->lineEdit_FWHMxy->setText(QString::number(gk.FWHM_xy));
            ui->lineEdit_FWHMz->setText(QString::number(gk.FWHM_z));
            ui->comboBox_FWHM_xy->setCurrentText(gk.unitFWHM_xy);
            ui->comboBox_FWHM_z->setCurrentText(gk.unitFWHM_z);
        }

    }
    if(!foundIt)
    {
        ui->lineEdit_FWHMxy->setText("");
        ui->lineEdit_FWHMz->setText("");
        ui->comboBox_FWHM_xy->setCurrentText("");
        ui->comboBox_FWHM_z->setCurrentText("");
    }

}

void sdmixer::SameKernelCheckBoxChanged(int state)
{
    if(ui->checkBox_sameKernel->isChecked())
    {
        ui->comboBox_ConvolutionChannel->setEnabled(false);
        ui->pushButton_saveChannel->setVisible(false);
    }
    else
    {
        ui->comboBox_ConvolutionChannel->setEnabled(true);
        ui->pushButton_saveChannel->setVisible(true);
    }
}

void sdmixer::loadKernel(int a)
{
qDebug() << "loadKernel";
    FilterFiles.clear();
    //qDeleteAll(ui->comboBox_ConvolutionChannel->set);
    ui->comboBox_ConvolutionChannel->clear();
    //qDebug() << ui->listWidget_FilterFiles->count();
    for(int i = 0; i < ui->listWidget_FilterFiles->count(); ++i)
    {
        QListWidgetItem* item = ui->listWidget_FilterFiles->item(i);
        qDebug()<< item->text();
        FilterFiles.push_back(item->text());
        QFileInfo fi(item->text());

        ui->comboBox_ConvolutionChannel->addItem(fi.baseName());
    }


    std::vector<gaussian_kernel>::iterator it=vec_kernel.begin();
    for( ; it != vec_kernel.end();)
    {
        bool stillExists=false;
        gaussian_kernel gk = *it;
        for(int i = 0; i < ui->comboBox_ConvolutionChannel->count(); ++i)
        {
            QString item = ui->comboBox_ConvolutionChannel->itemText(i);
            //qDebug()<< item;
            QFileInfo fi(item);
            QString baseName = fi.baseName();
            if(!gk.filterName.compare(baseName))
            {
                stillExists=true;
            }
        }
        if(!stillExists)
            vec_kernel.erase(it);
        else
            ++it;
    }

}

void sdmixer::UpdateShortPositionComboBox(QString str)
{
    //if(str.compare(ui->comboBox_CameraOrientation->currentText()))
    {
        ui->comboBox_ShortChannelPosition->clear();
        if(!str.compare("Left-Right"))
        {
            ui->comboBox_ShortChannelPosition->addItem("Left");
            ui->comboBox_ShortChannelPosition->addItem("Right");
            ui->comboBox_ShortChannelPosition->setCurrentIndex(0);
        }
        else
        {
            ui->comboBox_ShortChannelPosition->addItem("Top");
            ui->comboBox_ShortChannelPosition->addItem("Bottom");
            ui->comboBox_ShortChannelPosition->setCurrentIndex(0);
        }
    }
}

void sdmixer::CameraOrientationChanged(QString str)
{
    UpdateShortPositionComboBox(str);
}
void sdmixer::filterInfo(QString str)
{
    maxIntensityLong =  QString(ui->lineEdit_maxIntLong->text()).toDouble();
    maxIntensityShort = QString(ui->lineEdit_maxIntShort->text()).toDouble();
    precision = QString(ui->lineEdit_precision->text()).toDouble();
    FilterOrientation = ui->comboBox_FilterOrientation->currentText();
    std::stringstream ss;
    if (!FilterOrientation.compare("x=Short, y=Long"))
        ss <<"Filter files must be " << round(maxIntensityShort*precision) << " x " <<
         round(maxIntensityLong*precision) << " px";
    else
        ss <<"Filter files must be " << round(maxIntensityLong*precision) << " x " <<
         round(maxIntensityShort*precision) << " px";

    ui->label_FilterFileInstruction->setText(QString::fromStdString(ss.str()));
}


void sdmixer::InputFileClicked(QListWidgetItem *item)
{
    if(!item)
        return;

    PairFinder pf(this, item->text());
    ui->lineEdit_FileInfoDimension->setText(QString::number(pf.getDimensions()));
    ui->lineEdit_FileInfo_xmin->setText(QString::number(pf.getMinMaxValues().min_x*1e9));
    ui->lineEdit_FileInfo_xmax->setText(QString::number(pf.getMinMaxValues().max_x*1e9));
    ui->lineEdit_FileInfo_ymin->setText(QString::number(pf.getMinMaxValues().min_y*1e9));
    ui->lineEdit_FileInfo_ymax->setText(QString::number(pf.getMinMaxValues().max_y*1e9));
    ui->lineEdit_FileInfo_zmin->setText(QString::number(pf.getMinMaxValues().min_z*1e9));
    ui->lineEdit_FileInfo_zmax->setText(QString::number(pf.getMinMaxValues().max_z*1e9));
    getColsAndRows(item->text());
    ui->lineEdit_FileInfoCols->setText(QString::number(rawDataCols));
    ui->lineEdit_FileInfoRows->setText(QString::number(rawDataRows));

}

void sdmixer::getColsAndRows(QString file)
{
    rawDataCols=0;
    rawDataRows=0;
    std::ifstream ifs(file.toStdString());
    std::string firstLine, secondLine;
    getline (ifs, firstLine);
    getline (ifs, secondLine);
    ifs.close();
    //determine column number from second line
    std::stringstream countColsStream(secondLine);
    double dd;
    while (countColsStream >> dd)
    {
        ++rawDataCols;
    }
    // This is the fastest way to count lines in file
    // from the UNIX tool wget

    const char *fname = file.toStdString().c_str();

    static const int BUFFER_SIZE = 16*1024;
    int fd = open(fname, O_RDONLY);
    if(fd == -1)
        return;
    char buf[BUFFER_SIZE + 1];
    while(size_t bytes_read = read(fd, buf, BUFFER_SIZE))
    {
        if(bytes_read == (size_t)-1)
            return;//handle_error("read failed");
        if (!bytes_read)
            break;
        for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p)
            ++rawDataRows;
    }

    // Lines counted

}

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

    output_directory = ui->lineEdit_outputDirectory->text();

    pixelSizeNM = QString(ui->lineEdit_pixelSizeNM->text()).toInt();

    offset[0] = ui->lineEdit_xOffset->text().toDouble();
    offset[1] = ui->lineEdit_yOffset->text().toDouble();
    offset[2] = ui->lineEdit_zOffset->text().toDouble();
    epsilon[0] = ui->lineEdit_xEpsilon->text().toDouble();
    epsilon[1] = ui->lineEdit_yEpsilon->text().toDouble();
    epsilon[2] = ui->lineEdit_zEpsilon->text().toDouble();



    fishingSettings.run = ui->checkBox_runFishing->isChecked() ;
    fishingSettings.increment = QString(ui->lineEdit_FishingIncrement->text()).toDouble();
    fishingSettings.range = QString(ui->lineEdit_FishingRange->text()).toInt();
    fishingSettings.subset = QString(ui->lineEdit_FishingSubset->text()).toInt();
    fishingSettings.unit = ui->comboBox_fishingIncrement->currentText();


    offsetUnits.xOffset = ui->comboBox_xOffset->currentText();
    offsetUnits.yOffset = ui->comboBox_yOffset->currentText();
    offsetUnits.zOffset = ui->comboBox_zOffset->currentText();

    offsetUnits.xEpsilon = ui->comboBox_xEpsilon->currentText();
    offsetUnits.yEpsilon = ui->comboBox_yEpsilon->currentText();
    offsetUnits.zEpsilon = ui->comboBox_zEpsilon->currentText();
    CameraOrientation = ui->comboBox_CameraOrientation->currentText();
    ShortChannelPosition = ui->comboBox_ShortChannelPosition->currentText();

    // Filter
    FilterFiles.clear();
    for(int i = 0; i < ui->listWidget_FilterFiles->count(); ++i)
    {
        QListWidgetItem* item = ui->listWidget_FilterFiles->item(i);
        FilterFiles.push_back(item->text());
    }
    maxIntensityLong =  QString(ui->lineEdit_maxIntLong->text()).toDouble();
    maxIntensityShort = QString(ui->lineEdit_maxIntShort->text()).toDouble();
    precision = QString(ui->lineEdit_precision->text()).toDouble();
    FilterOrientation = ui->comboBox_FilterOrientation->currentText();


    // Reconstructor
    xyBinning = QString(ui->lineEdit_xyBinning->text()).toDouble();
    zBinning = QString(ui->lineEdit_zBinning->text()).toDouble();
    runConvolution = ui->checkBox_runConvolution->isChecked();
    //Convolution Kernel
    oneKernelForAllChannels = ui->checkBox_sameKernel->isChecked();
    if(oneKernelForAllChannels)
    {
        globalKernel.filterName = "global";
        globalKernel.FWHM_xy = QString(ui->lineEdit_FWHMxy->text()).toDouble();
        globalKernel.FWHM_z = QString(ui->lineEdit_FWHMz->text()).toDouble();
        globalKernel.unitFWHM_xy = ui->comboBox_FWHM_xy->currentText();
        globalKernel.unitFWHM_z = ui->comboBox_FWHM_z->currentText();
    }
    //loadKernel(ui->comboBox_ConvolutionChannel->currentIndex());

    nonLinearHistogramEqual = ui->checkBox_HistEq->isChecked();
    histeqCoefficient = QString(ui->lineEdit_CorrectionCoefficient->text()).toDouble();
    Threshold = QString(ui->lineEdit_Threshold->text()).toDouble();
    sqrtCummulation = ui->checkBox_sqrtCum->isChecked();



    LZWCompression = ui->checkBox_LZWCompression->isChecked();
    ResliceZ = ui->checkBox_scliceZ->isChecked();
    startRescliceZ = QString(ui->lineEdit_startSliceZ->text()).toInt();
    endRescliceZ = QString(ui->lineEdit_endSliceZ->text()).toInt();


    return true;

}
void sdmixer::setSettingsToUI(Settings s){

    // Session & Pairfinder
    std::vector<QString> vec = s.getInputFiles();
    for( auto i : vec)
    {
        insertItem(i, listWidgetInputFiles);
    }
    ui->lineEdit_outputDirectory->setText(s.getOutputDirectory());
    ui->lineEdit_pixelSizeNM->setText(QString::number(s.getPixelSizeNM()));
    ui->checkBox_runPairFinder->setChecked(s.getRunPairFinder());
    ui->checkBox_runFilter->setChecked(s.getRunFilter());
    ui->checkBox_runReconstructor->setChecked(s.getRunReconstructor());
    ui->checkBox_force2D->setChecked(s.getForce2D());


    ui->lineEdit_xOffset->setText(QString::number(s.getOffset(0)));
    ui->comboBox_xOffset->setCurrentText(s.getOffsetUnits().xOffset);
    ui->lineEdit_yOffset->setText(QString::number(s.getOffset(1)));
    ui->comboBox_yOffset->setCurrentText(s.getOffsetUnits().yOffset);
    ui->lineEdit_zOffset->setText(QString::number(s.getOffset(2)));
    ui->comboBox_zOffset->setCurrentText(s.getOffsetUnits().zOffset);

    ui->lineEdit_xEpsilon->setText(QString::number(s.getEpsilon(0)));
    ui->comboBox_xEpsilon->setCurrentText(s.getOffsetUnits().xEpsilon);
    ui->lineEdit_yEpsilon->setText(QString::number(s.getEpsilon(1)));
    ui->comboBox_yEpsilon->setCurrentText(s.getOffsetUnits().yEpsilon);
    ui->lineEdit_zEpsilon->setText(QString::number(s.getEpsilon(2)));
    ui->comboBox_zEpsilon->setCurrentText(s.getOffsetUnits().zEpsilon);

    fishing f =s.getFishing();

    ui->checkBox_runFishing->setChecked(f.run);
    ui->lineEdit_FishingIncrement->setText(QString::number(f.increment));
    ui->lineEdit_FishingRange->setText(QString::number(f.range));
    ui->lineEdit_FishingSubset->setText(QString::number(f.subset));

    ui->comboBox_fishingIncrement->setCurrentText(f.unit);
    ui->comboBox_CameraOrientation->setCurrentText(s.getCameraOrientation());
    ui->comboBox_ShortChannelPosition->clear();
    if(!s.getCameraOrientation().compare("Left-Right"))
    {
        ui->comboBox_ShortChannelPosition->addItem("Left");
        ui->comboBox_ShortChannelPosition->addItem("Right");
        ui->comboBox_ShortChannelPosition->setCurrentText(s.getShortChannelPosition());
    }else
    {
        ui->comboBox_ShortChannelPosition->addItem("Top");
        ui->comboBox_ShortChannelPosition->addItem("Bottom");
        ui->comboBox_ShortChannelPosition->setCurrentText(s.getShortChannelPosition());
    }

    // Filter
    std::vector<QString> v = s.getFilterFiles();
    for( auto i : v)
    {
        insertItem(i, listWidgetFilters);
    }
    ui->lineEdit_maxIntLong->setText(QString::number(s.getMaxIntLong()));
    ui->lineEdit_maxIntShort->setText(QString::number(s.getMaxIntShort()));
    ui->lineEdit_precision->setText(QString::number(s.getPrecision()));
    ui->comboBox_FilterOrientation->setCurrentText(s.getFilterOrientation());
    ui->comboBox_FilterOrientation->setCurrentText(s.getFilterOrientation());

    //Reconstructor
    ui->lineEdit_xyBinning->setText(QString::number(s.getXYbinning()));
    ui->lineEdit_zBinning->setText(QString::number(s.getZbinning()));
    ui->checkBox_runConvolution->setChecked(s.getRunConvolution());

    ui->checkBox_sameKernel->setChecked(s.getOneKernelForAllChannels());

    vec_kernel.clear();
    vec_kernel=s.getConvolutionKernel();
    gaussian_kernel gk;
    if(!vec_kernel.empty() && !s.getOneKernelForAllChannels())
        gk = vec_kernel[0];
    else
        gk = s.getGlobalKernel();

    ui->comboBox_ConvolutionChannel->setCurrentText(gk.filterName);
    ui->lineEdit_FWHMxy->setText(QString::number(gk.FWHM_xy));
    ui->lineEdit_FWHMz->setText(QString::number(gk.FWHM_z));
    ui->comboBox_FWHM_xy->setCurrentText(gk.unitFWHM_xy);
    ui->comboBox_FWHM_z->setCurrentText(gk.unitFWHM_z);

    ui->checkBox_HistEq->setChecked(s.getNonLinearHistEq());
    ui->lineEdit_CorrectionCoefficient->
            setText(QString::number(s.getCorrectionCoefficient()));
    ui->lineEdit_Threshold->setText(QString::number(s.getThreshold()));
    ui->checkBox_LZWCompression->setChecked(s.getLZWCompression());
    ui->checkBox_scliceZ->setChecked(s.getResliceZ());
    ui->lineEdit_startSliceZ->setText(QString::number(s.getStartSliceZ()));
    ui->lineEdit_endSliceZ->setText(QString::number(s.getEndSliceZ()));

}

void sdmixer::setStartDemixingButtonEnabled(bool val)
{
    ui->startDemixing->setEnabled(val);
}
void sdmixer::setCancelRunButtonEnabled(bool val)
{
    ui->pushButton_CancelRun->setVisible(val);
}

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
        QFileInfo fi(fileName);
        if(!fi.completeSuffix().compare("png"))
            insertItem(fileName, listWidgetFilters);
        else
            insertItem(fileName, listWidgetInputFiles);

    }
}
void sdmixer::insertItem(QString filename, QListWidget *list){
    QListWidgetItem *item = new QListWidgetItem(QIcon(filename), filename, list);
    list->setCurrentItem(item);
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

void sdmixer::threadReady()
{
    qDebug() << "thread Ready";
    nextStage();
}

void sdmixer::nextStage()
{

    ++current_stage;
    qDebug() << "current_stage:" << current_stage;
    if(current_stage == 1)
    {
        if ( runPairFinder )
        {
            QThread *pairfinder_thread = new QThread;
            PairFinder *pf;
            pf = new PairFinder(this, InputFiles[current_file]);
            //PairFinder p(this, InputFiles[current_file]);

            pf->moveToThread(pairfinder_thread);
            connect(pairfinder_thread, SIGNAL(started()), pf, SLOT(doWork()));
            connect(pf, SIGNAL(finished()), pairfinder_thread, SLOT(quit()));
            connect(pf, SIGNAL(finished()), pf, SLOT(deleteLater()));
            connect(pairfinder_thread, SIGNAL(finished()), pairfinder_thread, SLOT(deleteLater()));
            connect(pf, SIGNAL(finished()), this, SLOT(threadReady()));
            pairfinder_thread->start();
        }
        return;

    }
    if(current_stage == 2)
    {
        if ( runFilter )
        {
            QThread *filter_thread = new QThread;
            Filter *f;
            if( runPairFinder )
                f = new Filter(this, pf_output);
            else
                f = new Filter(this, InputFiles[current_file]);

            f->moveToThread(filter_thread);
            connect(filter_thread, SIGNAL(started()), f, SLOT(doWork()));
            connect(f, SIGNAL(finished()), filter_thread, SLOT(quit()));
            connect(f, SIGNAL(finished()), f, SLOT(deleteLater()));
            connect(filter_thread, SIGNAL(finished()), filter_thread, SLOT(deleteLater()));
            connect(f, SIGNAL(finished()),this, SLOT(threadReady()));
            filter_thread->start();
        }
        return;
    }
    if (current_stage == 3)
    {
        /*for(auto i:pf_output)
            qDebug() << i.filter;
        if( runReconstructor )
        {
            QThread *reconstructor_thread = new QThread;
            Reconstructor *r;

            if(runPairFinder && !runFilter) // Take all
                r = new Reconstructor(this, pf_output, 0);
            if(runFilter)
            {
                if(current_filter <= filter_max)
                    r = new Reconstructor(this, pf_output, current_filter);
            }
            if(!runPairFinder && !runFilter) // process xyz to simple images
                r = new Reconstructor(this, InputFiles[current_file]);

            r->moveToThread(reconstructor_thread);
            connect(reconstructor_thread, SIGNAL(started()), r, SLOT(doWork()));
            connect(r, SIGNAL(finished()), reconstructor_thread, SLOT(quit()));
            connect(r, SIGNAL(finished()), r, SLOT(deleteLater()));
            connect(reconstructor_thread, SIGNAL(finished()), reconstructor_thread, SLOT(deleteLater()));
            connect(r, SIGNAL(finished()),this, SLOT(threadReady()));
            reconstructor_thread->start();
        }*/
        return;

    }
    if (current_stage >= 4)
    {
        current_stage=0;

        current_file++;
        if( current_file >= InputFiles.size() )
        {
            qDebug() << current_file;
            current_stage=0;
            current_file = 0;
            sdmixer::log(this, "demixing finished!");
            setStartDemixingButtonEnabled(true);
            setCancelRunButtonEnabled(false);
        }
        else
        {
            current_stage = 0;
            qDebug() << "next stage";
            nextStage();
        }

        return;

    }
}


void sdmixer::on_startDemixing_clicked()
{
    getSettingsFromUI();

    if(InputFiles.empty())
    {
        QString err = error_msg("no input files selected");
        writeToConsole(err);
        return;
    }
    output_directory = ui->lineEdit_outputDirectory->text();
    ui->pushButton_CancelRun->setVisible(true);
    setStartDemixingButtonEnabled(false);

    nextStage();
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
    //ui->listWidget_FilterFiles->setCurrentRow(0);
    loadKernel(0);
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

void sdmixer::on_pushButton_CancelRun_clicked()
{
    ui->pushButton_CancelRun->setVisible(false);
    setStartDemixingButtonEnabled(true);
}

void sdmixer::on_pushButton_selectOutputDirectory_clicked()
{
    QString path = QFileDialog::getExistingDirectory(this, "Select Output Directory", QString());
    output_directory = path;
    ui->lineEdit_outputDirectory->setText(path);

}


void sdmixer::on_pushButton_saveChannel_clicked()
{

    if(! vec_kernel.empty())
    {
        std::vector<gaussian_kernel>::iterator it = vec_kernel.begin();
        for( ; it != vec_kernel.end();)
        {
            if(!it->filterName.compare(ui->comboBox_ConvolutionChannel->currentText()))
            {
                gaussian_kernel gk = *it;
                qDebug() <<" found it " << gk.filterName;
                vec_kernel.erase(it);
                gk.filterName=ui->comboBox_ConvolutionChannel->currentText();
                gk.FWHM_xy=ui->lineEdit_FWHMxy->text().toDouble();
                gk.FWHM_z=ui->lineEdit_FWHMz->text().toDouble();
                gk.unitFWHM_xy=ui->comboBox_FWHM_xy->currentText();
                gk.unitFWHM_z=ui->comboBox_FWHM_z->currentText();
                qDebug() << "not empty saved " << gk.filterName << " " << gk.FWHM_xy << " and " << gk.FWHM_z;
                vec_kernel.push_back(gk);
                //pushBackKernel(gk);

                qDebug() << "test";
                return;
            }
            else
            {
                it++;
            }
        }
    }

        qDebug() << ui->comboBox_ConvolutionChannel->currentText();
        gaussian_kernel gk2;
        gk2.filterName=ui->comboBox_ConvolutionChannel->currentText();
        gk2.FWHM_xy=ui->lineEdit_FWHMxy->text().toDouble();
        gk2.FWHM_z=ui->lineEdit_FWHMz->text().toDouble();
        gk2.unitFWHM_xy=ui->comboBox_FWHM_xy->currentText();
        gk2.unitFWHM_z=ui->comboBox_FWHM_z->currentText();
        qDebug() << "empty saved " << gk2.filterName << " " << gk2.FWHM_xy << " and " << gk2.FWHM_z << " and " << gk2.unitFWHM_xy << " and " << gk2.unitFWHM_z;


        vec_kernel.push_back(gk2);
        //pushBackKernel(gk2);

    qDebug() << "test2";
}
