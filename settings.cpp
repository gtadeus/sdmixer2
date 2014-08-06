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
    setInputFiles(s->getInputFiles());

    setPixelSizeNM(s->getPixelSize());

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

    setMaxIntLong(s->getMaxIntLong());
    setMaxIntShort(s->getMaxIntShort());
    setPrecision(s->getPrecision());

    setFilterFiles(s->getFilterFiles());

    setXYBinning(s->getReconstructor_xyBinning());
    setZBinning(s->getReconstructor_zBinning());

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
        qDebug()<<i;
        QDomElement child = appendChildNode(qdOffset, dimNames[i], offset[i]);
        child.setAttribute("unit", "px");
    }

    QDomElement qdEpsilon = createField("Epsilon");
    pairfinder.appendChild(qdEpsilon);

    for( int i = 0; i < max_dims; ++i)
    {
       QDomElement child = appendChildNode(qdEpsilon, dimNames[i], epsilon[i]);
       child.setAttribute("unit", "px");
    }

    QDomElement qdFishing = createField("Auto-OffsetSettings");
    pairfinder.appendChild(qdFishing);
    appendChildNode(qdFishing, "run", fishing.run);
    appendChildNode(qdFishing, "increment", fishing.increment);
    appendChildNode(qdFishing, "range", fishing.range);
    appendChildNode(qdFishing, "subset", fishing.subset);

    QDomElement filter = createField("FilterSettings");
    settings.appendChild(filter);
    appendChildNode(filter, "maxIntShort", maxIntensityShort);
    appendChildNode(filter, "maxIntLong", maxIntensityLong);
    appendChildNode(filter, "precision", precision);

    QDomElement reconstructor = createField("ReconstructorSettings");
    settings.appendChild(reconstructor);
    appendChildNode(reconstructor, "xyBinning", xyBinning);
    appendChildNode(reconstructor, "zBinning", zBinning);
    appendChildNode(reconstructor, "nonLinearHistEqual", false);
    appendChildNode(reconstructor, "runConvolution", runConvolution);
    appendChildNode(reconstructor, "FWHM_xy", 0);
    appendChildNode(reconstructor, "FWHM_z", 0);

    QDomElement expertSettings = createField("ExpertSettings");
    settings.appendChild(expertSettings);
    appendChildNode(expertSettings, "correctionCoefficient", 0.8);
    appendChildNode(expertSettings, "Threshold", 0);
    appendChildNode(expertSettings, "sqrtCummulation", true);

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
        if( e.tagName() == "field" )
        {
            if( e.attribute("name") == "GeneralSettings")
            {
                for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
                {
                    QDomElement f = m.toElement();
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
                            if(ee.attribute("name") == "x")
                                offset[0]=ee.attribute("number").toDouble();
                            if(ee.attribute("name") == "y")
                                offset[1]=ee.attribute("number").toDouble();
                            if(ee.attribute("name") == "z")
                                offset[2]=ee.attribute("number").toDouble();
                        }
                    }
                    if (f.attribute("name") == "Epsilon")
                    {
                        for(QDomNode mm = f.firstChild(); !mm.isNull(); mm = mm.nextSibling())
                        {
                            QDomElement ee = mm.toElement();
                            if(ee.attribute("name") == "x")
                                epsilon[0]=ee.attribute("number").toDouble();
                            if(ee.attribute("name") == "y")
                                epsilon[1]=ee.attribute("number").toDouble();
                            if(ee.attribute("name") == "z")
                                epsilon[2]=ee.attribute("number").toDouble();
                        }

                    }
                    if (f.attribute("name") == "Auto-OffsetSettings")
                    {
                        for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
                        {
                            QDomElement f = m.toElement();
                            if (f.attribute("name") == "run")
                            {
                                fishing.run = f.attribute("number").toInt();
                            }
                            if (f.attribute("name") == "increment")
                            {
                                fishing.increment = f.attribute("number").toInt();
                            }
                            if (f.attribute("name") == "range")
                            {
                                fishing.range = f.attribute("number").toInt();
                            }
                            if (f.attribute("name") == "subset")
                            {
                                fishing.subset = f.attribute("number").toInt();
                            }
                        }
                    }
                }
            }
            if( e.attribute("name") == "FilterSettings")
            {
                for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
                {
                    QDomElement f = m.toElement();
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

                    }
                    if (f.attribute("name") == "runConvolution")
                    {
                        runConvolution = f.attribute("number").toInt();
                    }
                    if (f.attribute("name") == "FWHM_xy")
                    {

                    }
                    if (f.attribute("name") == "FWHM_z")
                    {
                    }
                }
            }
            if( e.attribute("name") == "ExpertSettings")
            {
                for(QDomNode m = e.firstChild(); !m.isNull(); m = m.nextSibling())
                {
                    QDomElement f = m.toElement();
                    if (f.attribute("name") == "correctionCoefficient")
                    {
                    }
                    if (f.attribute("name") == "Threshold")
                    {
                    }
                    if (f.attribute("name") == "sqrtCummulation")
                    {
                    }
                }
            }
        }
    }
}
