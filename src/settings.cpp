#include "settings.h"

void Settings::setInputFiles(std::vector<QString> str)
{
    if(str.empty())
    {
        qDebug("settings: no input files");
    }
    else
    {
        this->InputFiles=str;
    }
}

Settings::Settings(sdmixer *s)
{
    // Session & Pairfinder
    setInputFiles(s->getInputFiles());
    setPixelSizeNM(s->getPixelSize());
    output_directory = s->getOutputDirectory();

    setRunPairfinder(s->getRunPairfinder());
    setRunFilter(s->getRunFilter());
    setRunReconstructor(s->getRunReconstructor());

    setForce2D(s->getForce2D());

    for(int i = 0; i< max_dims; ++i)
    {
        setOffset(i, s->getOffset(i));
        setEpsilon(i, s->getEpsilon(i));
    }

    setFishing(s->getFishing());
    CameraOrientation = s->getCameraOrientation();
    ShortChannelPosition = s->getShortChannelPosition();
    offsetUnits = s->getOffsetUnits();

    runGrouping = s->getRunGrouping();
    groupingRadius = s->getGroupingRadius();
    groupingUnits = s->getGroupingUnits();


    // Filter
    setFilterFiles(s->getFilterFiles());
    setMaxIntLong(s->getMaxIntLong());
    setMaxIntShort(s->getMaxIntShort());
    setPrecision(s->getPrecision());
    FilterOrientation = s->getFilterOrientation();
    plotIntensitySpace = s->getPlotIntensitySpace();

    //Reconstructor
    setXYBinning(s->getReconstructor_xyBinning());
    setZBinning(s->getReconstructor_zBinning());
    runConvolution = s->getRunConvolution();
    oneKernelForAllChannels = s->getOneKernelForAllChannels();
    global_kernel = s->getGlobalKernel();
    vec_kernel = s->getConvolutionKernel();
    nonLinearHistogramEqual = s->getNonLinearHistEq();
    histeqCoefficient = s->getHisteqCoefficient();
    Threshold = s->getThreshold();
    sqrtCummulation = s->getSqrtCummulation();

    LZWCompression = s->getLZWCompression();
    ResliceZ = s->getResliceZ();
    startRescliceZ = s->getStartRescliceZ();
    endRescliceZ = s->getEndRescliceZ();

    performNNStatistic = s->getPerformNNStatistic();

}

void Settings::writeSettingsToFile(QString fileName)
{
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
            return;

    QTextStream out(&file);
    out << settingsFile.toString();

}
Settings::Settings(QString file){
    loadFromFile(file);
}
void Settings::initXML(){
    QDomProcessingInstruction header = settingsFile.createProcessingInstruction("xml", "version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ");
    settingsFile.appendChild(header);

    QDomElement settings = settingsFile.createElement("ConfigFile");
    settings.setAttribute("name", "sdmixer");
    settings.setAttribute("version", "2.01");
    settingsFile.appendChild(settings);

    if ( !InputFiles.empty() )
    {
        QDomElement qdInputFiles = settingsFile.createElement( "InputFiles" );
        settings.appendChild( qdInputFiles );
        for( auto i : InputFiles)
        {
            QDomElement qdInputFileEntry = settingsFile.createElement( "InputFile" );
            qdInputFileEntry.setAttribute( "path", i );
            qdInputFiles.appendChild( qdInputFileEntry );
        }
    }

    QDomElement general = settingsFile.createElement("field");
    general.setAttribute("name", "GeneralSettings");

    appendChildNode(general, "OutputDirectory", output_directory);
    appendChildNode(general, "pixelSizeNM", pixelSizeNM);
    appendChildNode(general, "runPairFinder", runPairFinder);
    appendChildNode(general, "runFilter", runFilter);
    appendChildNode(general, "runReconstructor", runReconstructor);
    appendChildNode(general, "force2D", force2D);

    settings.appendChild(general);
    QDomElement pairfinder = createField("PairFinderSettings");
    settings.appendChild(pairfinder);


    QDomElement qdOffset = createField("Offset");
    pairfinder.appendChild(qdOffset);

    QString dimNames[max_dims] = {"x", "y", "z"};

    for( int i = 0; i < max_dims; ++i)
    {
        QDomElement child = appendChildNode(qdOffset, dimNames[i], offset[i]);
        child.setAttribute("unit", offsetUnits.getOffset(i));
    }

    QDomElement qdEpsilon = createField("Epsilon");
    pairfinder.appendChild(qdEpsilon);

    for( int i = 0; i < max_dims; ++i)
    {
       QDomElement child = appendChildNode(qdEpsilon, dimNames[i], epsilon[i]);
       child.setAttribute("unit", offsetUnits.getEpsilon(i));
    }

    appendChildNode(pairfinder, "CameraOrientation", CameraOrientation);
    appendChildNode(pairfinder, "ShortChannelPosition", ShortChannelPosition);

    QDomElement qdFishing = createField("Auto-OffsetSettings");
    pairfinder.appendChild(qdFishing);
    appendChildNode(qdFishing, "run", fishing.run);
    appendChildNode(qdFishing, "increment", fishing.increment);
    appendChildNode(qdFishing, "unit", fishing.unit);
    appendChildNode(qdFishing, "range", fishing.range);
    appendChildNode(qdFishing, "subset", fishing.subset);

    appendChildNode(pairfinder, "runGrouping", runGrouping);
    appendChildNode(pairfinder, "groupingRadius", groupingRadius);
    appendChildNode(pairfinder, "groupingUnits", groupingUnits);


    QDomElement filter = createField("FilterSettings");
    settings.appendChild(filter);
    if ( !FilterFiles.empty() )
    {
        QDomElement qdFilterFiles = settingsFile.createElement( "FilterFiles" );
        filter.appendChild( qdFilterFiles );
        for( auto i : FilterFiles)
        {
            QDomElement qdFilterFileEntry = settingsFile.createElement( "FilterFile" );
            qdFilterFileEntry.setAttribute( "path", i );
            qdFilterFiles.appendChild( qdFilterFileEntry );
        }
    }
    appendChildNode(filter, "maxIntShort", maxIntensityShort);
    appendChildNode(filter, "maxIntLong", maxIntensityLong);
    appendChildNode(filter, "precision", precision);
    appendChildNode(filter, "FilterOrientation", FilterOrientation);
    appendChildNode(filter, "plotIntensitySpace", plotIntensitySpace);

    QDomElement reconstructor = createField("ReconstructorSettings");
    settings.appendChild(reconstructor);
    appendChildNode(reconstructor, "xyBinning", xyBinning);
    appendChildNode(reconstructor, "zBinning", zBinning);
    appendChildNode(reconstructor, "nonLinearHistEqual", nonLinearHistogramEqual);
    appendChildNode(reconstructor, "runConvolution", runConvolution);
    appendChildNode(reconstructor, "oneKernelForAllChannels", oneKernelForAllChannels);
    appendChildNode(reconstructor, "FWHM_xy", global_kernel.FWHM_xy);
    appendChildNode(reconstructor, "FWHM_z", global_kernel.FWHM_z);
    appendChildNode(reconstructor, "unitFWHM_xy", global_kernel.unitFWHM_xy);
    appendChildNode(reconstructor, "unitFWHM_z", global_kernel.unitFWHM_z);


    appendChildNode(reconstructor, "correctionCoefficient", histeqCoefficient);
    appendChildNode(reconstructor, "Threshold", Threshold);
    appendChildNode(reconstructor, "sqrtCummulation", sqrtCummulation);
    appendChildNode(reconstructor, "LZWCompression", LZWCompression);
    appendChildNode(reconstructor, "ResliceZ", ResliceZ);
    appendChildNode(reconstructor, "startRescliceZ", startRescliceZ);
    appendChildNode(reconstructor, "endRescliceZ", endRescliceZ);

    appendChildNode(reconstructor, "NearestNeighborStatistic", performNNStatistic);

    std::vector<sdmixer::gaussian_kernel>::iterator it;
    for( it = vec_kernel.begin(); it != vec_kernel.end(); ++it)
    {
        QDomElement qdConvKernel = createField("ConvolutionKernel");
        sdmixer::gaussian_kernel gk=*it;
        appendChildNode(qdConvKernel, "FilterFile", gk.filterName);
        appendChildNode(qdConvKernel, "FWHM_xy", gk.FWHM_xy);
        appendChildNode(qdConvKernel, "FWHM_z", gk.FWHM_z);
        appendChildNode(qdConvKernel, "FWHM_xy_unit", gk.unitFWHM_xy);
        appendChildNode(qdConvKernel, "FWHM_z_unit", gk.unitFWHM_z);
        reconstructor.appendChild(qdConvKernel);
    }


    //QDomElement expertSettings = createField("ExpertSettings");
    //settings.appendChild(expertSettings);

    settingsFile.appendChild(settings);
}
QDomElement Settings::createField(QString name){
    QDomElement qd = settingsFile.createElement("field");
    qd.setAttribute("name", name);
    return qd;
}

QDomElement Settings::appendChildNode(QDomElement parent, QString name, double val){
    QDomElement qd = settingsFile.createElement("value");
    qd.setAttribute("number", val);
    qd.setAttribute("name", name);
    parent.appendChild(qd);
    return qd;
}
QDomElement Settings::appendChildNode(QDomElement parent, QString name, QString str){
    QDomElement qd = settingsFile.createElement("value");
    qd.setAttribute("string", str);
    qd.setAttribute("name", name);
    parent.appendChild(qd);
    return qd;
}

void Settings::loadFromFile(QString file){
    QFile inFile( file );
    if( !inFile.open( QIODevice::ReadOnly | QIODevice::Text ) )
    {
      qDebug( "Failed to open file for reading." );
      //return 0;
    }

    if( !settingsFile.setContent( &inFile ) )
    {
      qDebug( "Failed to parse the file into a DOM tree." );
      inFile.close();

    }
    QDomElement element = settingsFile.documentElement();

    for(QDomNode n = element.firstChild(); !n.isNull(); n = n.nextSibling())
    {
        QDomElement e = n.toElement();
        if( e.tagName() == "InputFiles")
        {
            for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
            {
                QDomElement f = m.toElement();
                InputFiles.push_back(f.attribute("path"));
            }
        }
        if( e.tagName() == "field" )
        {
            if( e.attribute("name") == "GeneralSettings")
            {
                for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
                {
                    QDomElement f = m.toElement();
                    if (f.attribute("name") == "OutputDirectory")
                    {
                        output_directory = f.attribute("string");
                    }
                    if (f.attribute("name") == "pixelSizeNM")
                    {
                        pixelSizeNM = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "runPairFinder")
                    {
                        runPairFinder = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "runFilter")
                    {
                        runFilter = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "runReconstructor")
                    {
                        runReconstructor = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "force2D")
                    {
                        force2D = f.attribute("number").toInt();
                    }
                }
            }
            if( e.attribute("name") == "PairFinderSettings")
            {
                for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
                {
                    QDomElement f = m.toElement();
                    if (f.attribute("name") == "Offset")
                    {
                        for(QDomNode mm = f.firstChild(); !mm.isNull(); mm = mm.nextSibling())
                        {
                            QDomElement ee = mm.toElement();
                            if(ee.attribute("name") == "x"){
                                offset[0]=ee.attribute("number").toDouble();
                                offsetUnits.xOffset=ee.attribute("unit");
                            }
                            if(ee.attribute("name") == "y"){
                                offset[1]=ee.attribute("number").toDouble();
                                offsetUnits.yOffset=ee.attribute("unit");
                            }
                            if(ee.attribute("name") == "z"){
                                offset[2]=ee.attribute("number").toDouble();
                                offsetUnits.zOffset=ee.attribute("unit");
                            }
                        }
                    }
                    if (f.attribute("name") == "Epsilon")
                    {
                        for(QDomNode mm = f.firstChild(); !mm.isNull(); mm = mm.nextSibling())
                        {
                            QDomElement ee = mm.toElement();
                            if(ee.attribute("name") == "x"){
                                epsilon[0]=ee.attribute("number").toDouble();
                                offsetUnits.xEpsilon=ee.attribute("unit");
                            }
                            if(ee.attribute("name") == "y"){
                                epsilon[1]=ee.attribute("number").toDouble();
                                offsetUnits.yEpsilon=ee.attribute("unit");
                            }
                            if(ee.attribute("name") == "z"){
                                epsilon[2]=ee.attribute("number").toDouble();
                                offsetUnits.zEpsilon=ee.attribute("unit");
                            }
                        }

                    }
                    if (f.attribute("name") == "CameraOrientation")
                    {
                        CameraOrientation = f.attribute("string");
                    }
                    if (f.attribute("name") == "ShortChannelPosition")
                    {
                        ShortChannelPosition = f.attribute("string");
                    }
                    if (f.attribute("name") == "Auto-OffsetSettings")
                    {
                        for(QDomNode m = f.firstChild(); !m.isNull(); m = m.nextSibling())
                        {
                            QDomElement f = m.toElement();
                            if (f.attribute("name") == "run")
                            {
                                fishing.run = f.attribute("number").toInt();
                            }
                            if (f.attribute("name") == "increment")
                            {
                                fishing.increment = f.attribute("number").toDouble();
                            }
                            if (f.attribute("name") == "range")
                            {
                                fishing.range = f.attribute("number").toInt();
                            }
                            if (f.attribute("name") == "subset")
                            {
                                fishing.subset = f.attribute("number").toInt();
                            }
                            if (f.attribute("name") == "unit")
                            {
                                fishing.unit = f.attribute("string");
                            }
                        }
                    }
                    if (f.attribute("name") == "runGrouping")
                    {
                        runGrouping = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "groupingRadius")
                    {
                        groupingRadius = f.attribute("number").toDouble();
                    }
                    if (f.attribute("name") == "groupingUnits")
                    {
                        groupingUnits = f.attribute("string");
                    }
                }
            }
            if( e.attribute("name") == "FilterSettings")
            {
                for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
                {
                    QDomElement f = m.toElement();
                    if( f.tagName() == "FilterFiles")
                    {
                        for(QDomNode mm = m.firstChild(); !mm.isNull(); mm = mm.nextSibling())
                        {
                            QDomElement ff = mm.toElement();
                            FilterFiles.push_back(ff.attribute("path"));
                        }
                    }
                    if (f.attribute("name") == "maxIntShort")
                    {
                        maxIntensityShort = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "maxIntLong")
                    {
                        maxIntensityLong = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "precision")
                    {
                        precision = f.attribute("number").toDouble();
                    }
                    if (f.attribute("name") == "FilterOrientation")
                    {
                        FilterOrientation = f.attribute("string");
                    }
                    if (f.attribute("name") == "plotIntensitySpace")
                    {
                        plotIntensitySpace = f.attribute("number").toInt();
                    }
                }
            }
            if( e.attribute("name") == "ReconstructorSettings")
            {
                for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
                {
                    QDomElement f = m.toElement();
                    if (f.attribute("name") == "xyBinning")
                    {
                        xyBinning = f.attribute("number").toDouble();
                    }
                    if (f.attribute("name") == "zBinning")
                    {
                        zBinning = f.attribute("number").toDouble();
                    }
                    if (f.attribute("name") == "nonLinearHistEqual")
                    {
                        nonLinearHistogramEqual = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "runConvolution")
                    {
                        runConvolution = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "correctionCoefficient")
                    {
                        histeqCoefficient = f.attribute("number").toDouble();
                    }
                    if (f.attribute("name") == "Threshold")
                    {
                        Threshold = f.attribute("number").toDouble();
                    }
                    if (f.attribute("name") == "sqrtCummulation")
                    {
                        sqrtCummulation = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "LZWCompression")
                    {
                        LZWCompression = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "ResliceZ")
                    {
                        ResliceZ = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "startRescliceZ")
                    {
                        startRescliceZ = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "endRescliceZ")
                    {
                        endRescliceZ = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "oneKernelForAllChannels")
                    {
                        oneKernelForAllChannels = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "FWHM_xy")
                    {
                        global_kernel.FWHM_xy = f.attribute("number").toDouble();
                    }
                    if (f.attribute("name") == "FWHM_z")
                    {
                        global_kernel.FWHM_z = f.attribute("number").toDouble();
                    }
                    if (f.attribute("name") == "unitFWHM_xy")
                    {
                        global_kernel.unitFWHM_xy = f.attribute("string");
                    }
                    if (f.attribute("name") == "unitFWHM_z")
                    {
                        global_kernel.unitFWHM_z = f.attribute("string");
                    }
                    if(f.attribute("name") == "ConvolutionKernel")
                    {
                        sdmixer::gaussian_kernel gk;
                        QDomElement ff = f.toElement();
                        //qDebug() << "hier";
                        for(QDomNode mm = ff.firstChild(); !mm.isNull(); mm = mm.nextSibling())
                        {
                            //qDebug() << "hier";
                            QDomElement nn = mm.toElement();
                            if (nn.attribute("name") == "FilterFile")
                            {
                                gk.filterName = nn.attribute("string");
                                qDebug() << gk.filterName;
                            }
                            if (nn.attribute("name") == "FWHM_xy")
                            {
                                gk.FWHM_xy = nn.attribute("number").toDouble();
                            }
                            if (nn.attribute("name") == "FWHM_z")
                            {
                                gk.FWHM_z = nn.attribute("number").toDouble();
                            }
                            if (nn.attribute("name") == "FWHM_xy_unit")
                            {
                                gk.unitFWHM_xy = nn.attribute("string");
                            }
                            if (nn.attribute("name") == "FWHM_z_unit")
                            {
                                gk.unitFWHM_z = nn.attribute("string");
                            }
                        }
                        vec_kernel.push_back(gk);
                    }
                    if (f.attribute("name") == "NearestNeighborStatistic")
                    {
                        performNNStatistic = f.attribute("number").toInt();
                    }
                }
            }
            if( e.attribute("name") == "ExpertSettings")
            {
                for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
                {
                    //QDomElement f = m.toElement();
                }
            }
        }
    }
}
