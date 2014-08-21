#include "sdmixer.h"
#include <QApplication>
#include <QtDebug>
#include <QDir>
#include <stdio.h>
#include <stdlib.h>

//#define SDMIXER_DEBUG_FILE "sdmixer_debug_log.txt"
QString debug_log_file;

void myMessageOutput(QtMsgType type, const QMessageLogContext &context, const QString &msg)
{
    FILE *file;

    file = fopen(debug_log_file.toLocal8Bit().data(), "a");
    QByteArray localMsg = msg.toLocal8Bit();
    switch (type) {
    case QtDebugMsg:
        fprintf(file, "%s\n", localMsg.constData());
        fclose(file);
        //fprintf(file, "Debug: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        //fclose(file);
        break;
    case QtWarningMsg:
        fprintf(file, "Warning: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        fclose(file);
        break;
    case QtCriticalMsg:
        fprintf(file, "Critical: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        fclose(file);
        break;
    case QtFatalMsg:
        fprintf(file, "Fatal: %s\n", localMsg.constData());
        fclose(file);
        //fprintf(file, "Fatal: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        //fclose(file);
        abort();
    }
}
bool file_exists (char *filename)
{
    if (FILE *file = fopen(filename, "r"))
    {
        fclose(file);
        return true;
    }
    else
        return false;

}

int main(int argc, char *argv[])
{
    debug_log_file = QDir::homePath();
    debug_log_file.append("/sdmixer_debug_log.txt");

    if(file_exists(debug_log_file.toLocal8Bit().data()))
        remove(debug_log_file.toLocal8Bit().data());

    qInstallMessageHandler(myMessageOutput);
    QApplication a(argc, argv);
    sdmixer w;
    w.show();

    return a.exec();
}
