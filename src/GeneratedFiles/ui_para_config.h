/********************************************************************************
** Form generated from reading UI file 'para_config.ui'
**
** Created by: Qt User Interface Compiler version 5.1.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PARA_CONFIG_H
#define UI_PARA_CONFIG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_para_config
{
public:
    QWidget *formLayoutWidget;
    QFormLayout *formLayout_2;
    QLabel *label;
    QDoubleSpinBox *data;
    QLabel *label_2;
    QDoubleSpinBox *rigid;
    QLabel *label_3;
    QDoubleSpinBox *smooth;
    QLabel *lable6;
    QDoubleSpinBox *tolerance;
    QPushButton *okButton;

    void setupUi(QWidget *para_config)
    {
        if (para_config->objectName().isEmpty())
            para_config->setObjectName(QStringLiteral("para_config"));
        para_config->resize(541, 756);
        formLayoutWidget = new QWidget(para_config);
        formLayoutWidget->setObjectName(QStringLiteral("formLayoutWidget"));
        formLayoutWidget->setGeometry(QRect(10, 20, 192, 101));
        formLayout_2 = new QFormLayout(formLayoutWidget);
        formLayout_2->setObjectName(QStringLiteral("formLayout_2"));
        formLayout_2->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
        formLayout_2->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(formLayoutWidget);
        label->setObjectName(QStringLiteral("label"));

        formLayout_2->setWidget(0, QFormLayout::LabelRole, label);

        data = new QDoubleSpinBox(formLayoutWidget);
        data->setObjectName(QStringLiteral("data"));
        data->setDecimals(4);
        data->setMaximum(10000);

        formLayout_2->setWidget(0, QFormLayout::FieldRole, data);

        label_2 = new QLabel(formLayoutWidget);
        label_2->setObjectName(QStringLiteral("label_2"));

        formLayout_2->setWidget(1, QFormLayout::LabelRole, label_2);

        rigid = new QDoubleSpinBox(formLayoutWidget);
        rigid->setObjectName(QStringLiteral("rigid"));
        rigid->setDecimals(4);
        rigid->setMaximum(10000);

        formLayout_2->setWidget(1, QFormLayout::FieldRole, rigid);

        label_3 = new QLabel(formLayoutWidget);
        label_3->setObjectName(QStringLiteral("label_3"));

        formLayout_2->setWidget(2, QFormLayout::LabelRole, label_3);

        smooth = new QDoubleSpinBox(formLayoutWidget);
        smooth->setObjectName(QStringLiteral("smooth"));
        smooth->setDecimals(4);
        smooth->setMaximum(10000);

        formLayout_2->setWidget(2, QFormLayout::FieldRole, smooth);

        lable6 = new QLabel(formLayoutWidget);
        lable6->setObjectName(QStringLiteral("lable6"));

        formLayout_2->setWidget(3, QFormLayout::LabelRole, lable6);

        tolerance = new QDoubleSpinBox(formLayoutWidget);
        tolerance->setObjectName(QStringLiteral("tolerance"));
        tolerance->setDecimals(4);

        formLayout_2->setWidget(3, QFormLayout::FieldRole, tolerance);

        okButton = new QPushButton(para_config);
        okButton->setObjectName(QStringLiteral("okButton"));
        okButton->setGeometry(QRect(60, 120, 75, 23));

        retranslateUi(para_config);

        QMetaObject::connectSlotsByName(para_config);
    } // setupUi

    void retranslateUi(QWidget *para_config)
    {
        para_config->setWindowTitle(QApplication::translate("para_config", "Config", 0));
        label->setText(QApplication::translate("para_config", "<html><head/><body><p><span style=\" font-size:11pt;\">data</span></p></body></html>", 0));
        label_2->setText(QApplication::translate("para_config", "<html><head/><body><p><span style=\" font-size:11pt;\">rigid</span></p></body></html>", 0));
        label_3->setText(QApplication::translate("para_config", "<html><head/><body><p><span style=\" font-size:11pt;\">smooth</span></p></body></html>", 0));
        lable6->setText(QApplication::translate("para_config", "<html><head/><body><p>tolerance </p></body></html>", 0));
        okButton->setText(QApplication::translate("para_config", "OK", 0));
    } // retranslateUi

};

namespace Ui {
    class para_config: public Ui_para_config {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PARA_CONFIG_H
