/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "rendersettingsdlg.h"
#include "ui_rendersettingsdlg.h"
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>

class BetterSpinBox : public QSpinBox {
public:
	BetterSpinBox(QWidget *parent) : QSpinBox(parent) {
		setFrame(false);
		setFocusPolicy(Qt::StrongFocus);
		setMinimum(-1);
		setMaximum(INT_MAX);
	}

	bool event(QEvent *event) {
		if(event->type() == QEvent::Wheel) {
			event->ignore();
			return false;
		}
		return QSpinBox::event(event);
	}
};

class BetterDoubleSpinBox : public QDoubleSpinBox {
public:
	BetterDoubleSpinBox(QWidget *parent) : QDoubleSpinBox(parent) {
		setDecimals(5);
		setMinimum(-std::numeric_limits<double>::max());
		setMaximum(std::numeric_limits<double>::max());
		setFrame(false);
		setFocusPolicy(Qt::StrongFocus);
	}

	void morphNumericString(char *s) const {
		char *p = strchr (s,'.');
		if (p == NULL) {
			strcat(s, ".0");
			return;
		}
		p = &(p[strlen(p)-1]);
		while ((p != s) && (*p == '0') && (*(p-1) != '.'))
			*p-- = '\0';
	}

	QString textFromValue(double value) const {
		char tmp[32];
		snprintf(tmp, sizeof(tmp), "%f", value);
		morphNumericString(tmp);
		return QString::fromAscii(tmp);
	}

	bool event(QEvent *event) {
		if(event->type() == QEvent::Wheel) {
			event->ignore();
			return false;
		}
		return QDoubleSpinBox::event(event);
	}
};

/* ====================== Some helper routines ====================== */

static void setComboBox(QComboBox *box, const std::string &pluginName) {
	for (int i=0; i<box->count(); ++i) {
		const QList<QVariant> &data = box->itemData(i).toList();
		if (data.at(2) == pluginName.c_str()) {
			box->setCurrentIndex(i);
			return;
		}
	}
	SLog(EError, "The plugin \"%s\" can't be controlled via the GUI! (update <tt>"
		"src/mtsgui/resources/docs.xml</tt> and recompile in case you want to be "
		"able to do this)", pluginName.c_str());
}

static std::string getPluginName(QComboBox *box) {
	return box->itemData(box->currentIndex()).toList().at(2).toString().toStdString();
}

/* ====================== RenderSettingsDialog impl ====================== */

RenderSettingsDialog::RenderSettingsDialog(QWidget *parent) :
		QDialog(parent, Qt::Sheet),
	ui(new Ui::RenderSettingsDialog), m_icNode(NULL), m_aiNode(NULL) {
	ui->setupUi(this);

	connect(ui->integratorBox, SIGNAL(highlighted(int)), SLOT(cbHighlighted(int)));
	connect(ui->integratorBox, SIGNAL(activated(int)), SLOT(update()));
	connect(ui->samplerBox, SIGNAL(highlighted(int)), SLOT(cbHighlighted(int)));
	connect(ui->samplerBox, SIGNAL(activated(int)), SLOT(update()));
	connect(ui->rFilterBox, SIGNAL(highlighted(int)), SLOT(cbHighlighted(int)));
	connect(ui->rFilterBox, SIGNAL(activated(int)), SLOT(update()));
	connect(ui->icBox, SIGNAL(pressed()), SLOT(chkBoxPressed()));
	connect(ui->aiBox, SIGNAL(pressed()), SLOT(chkBoxPressed()));
	connect(ui->icBox, SIGNAL(toggled(bool)), SLOT(update()));
	connect(ui->aiBox, SIGNAL(toggled(bool)), SLOT(update()));
	connect(ui->resolutionBox, SIGNAL(activated(int)), SLOT(refresh()));
	connect(ui->resolutionBox, SIGNAL(editTextChanged(const QString &)), SLOT(refresh()));

	QFile file(":/resources/docs.xml");
	if (!file.open(QIODevice::ReadOnly) || !m_document.setContent(&file))
		SLog(EError, "Unable to read the documentation file!");
	file.close();

	/* Populate the integrator, rec. filter & sampler combo box widgets */
	QDomElement docRoot = m_document.documentElement();
	for (QDomElement e = docRoot.firstChildElement("plugin"); !e.isNull();
		 e = e.nextSiblingElement("plugin")) {
		QString docString, name = e.attribute("name");
		if (!e.firstChildElement("descr").isNull()) {
			/* Create a HTML-based documentation string */
			QDomDocument helpDoc;
			QDomElement root = helpDoc.createElement("p");
			helpDoc.appendChild(root);

			for (QDomNode child = e.firstChildElement("descr").firstChild();
			   !child.isNull(); child = child.nextSibling())
			root.appendChild(helpDoc.importNode(child, true));
			docString = helpDoc.toString();
		}

		if (e.hasAttribute("show") && e.attribute("show") == "true") {
			QString type = e.attribute("type"),
					className = e.attribute("className"),
					readableName = e.attribute("readableName"),
					name = e.attribute("name");

			QList<QVariant> list;
			list.append(className);
			list.append(docString);
			list.append(name);

			if (type == "integrator")
				ui->integratorBox->addItem(readableName, list);
			else if (type == "sampler")
				ui->samplerBox->addItem(readableName, list);
			else if (type == "rfilter")
				ui->rFilterBox->addItem(readableName, list);
		}
		if (name == "irrcache")
			ui->icBox->setProperty("help", docString);
		else if (name == "adaptive")
			ui->aiBox->setProperty("help", docString);
	}

	m_model = new XMLTreeModel(docRoot, palette(), this);
	ui->treeView->setModel(m_model);
	ui->treeView->setAlternatingRowColors(true);
	ui->treeView->setUniformRowHeights(true);
	ui->treeView->setColumnWidth(0, 270);
	ui->treeView->setStyleSheet("QTreeView::item { padding-right: 8px; }");
	ui->treeView->setItemDelegate(new PropertyDelegate(this));
	connect(ui->treeView->selectionModel(), SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection)),
		SLOT(onTreeSelectionChange(const QItemSelection &, const QItemSelection &)));
	connect(m_model, SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(dataChanged()));
	m_integratorNode = m_model->registerClass("MIPathTracer", "Path tracer");
	m_samplerNode = m_model->registerClass("IndependentSampler", "Independent sampler");
	m_rFilterNode = m_model->registerClass("BoxFilter", "Box filter");
	QRegExp resRegExp("^[1-9]\\d{0,4}x[1-9]\\d{0,4}$");
	ui->resolutionBox->setValidator(new QRegExpValidator(resRegExp, this));
	QPalette pal = ui->helpViewer->palette();
	pal.setColor(QPalette::Text, pal.color(QPalette::Foreground));
	pal.setColor(QPalette::Base, pal.color(QPalette::Window));
	ui->helpViewer->setPalette(pal);
	ui->helpViewer->setHtml("Click on any setting for documentation");
}

void RenderSettingsDialog::setDocumentation(const QString &text) {
	m_currentDocumentation = text;
	bool hasErrors = false;
	QString comments;

	if (m_statusMessages.size() > 0) {
		comments = QString("<ul style=\"margin:0px;\">");
		ui->groupBox->setTitle(tr("Issues with the current configuration"));
		for (int i=0; i<m_statusMessages.size(); ++i) {
			const QString &message = m_statusMessages[i];
			bool isWarning = false, isError = false;

			if (message.contains("Warning"))
				isWarning = true;
			if (message.contains("Error"))
				isError = true;

			if (isError || isWarning)
				comments += QString("<li><em>%1</em></li>").arg(message);
			else
				comments += QString("<li>%1</li>").arg(message);
			hasErrors |= isError;
		}
		comments += "</ul>";
	} else {
		ui->groupBox->setTitle(tr("Documentation"));
	}

	#if defined(__OSX__)
		ui->helpViewer->setHtml(comments + "<div style='font-size:12pt'>" + m_currentDocumentation + "</div>");
	#else
		ui->helpViewer->setHtml(comments + "<div style='font-size:10pt'>" + m_currentDocumentation + "</div>");
	#endif
   	ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(!hasErrors);
}

void RenderSettingsDialog::dataChanged() {
	QStringList statusMessages = validateConfiguration();
	if (statusMessages != m_statusMessages) {
		m_statusMessages = statusMessages;
		setDocumentation(m_currentDocumentation);
	}
}

void RenderSettingsDialog::cbHighlighted(int index) {
	QComboBox *comboBox = static_cast<QComboBox *>(sender());
	setDocumentation(comboBox->itemData(index).toList().at(1).toString());
}

void RenderSettingsDialog::chkBoxPressed() {
	QCheckBox *checkBox = static_cast<QCheckBox *>(sender());
	setDocumentation(checkBox->property("help").toString());
}

void RenderSettingsDialog::onTreeSelectionChange(const QItemSelection &selected, const QItemSelection &deselected) {
	QModelIndexList indexList = selected.indexes();
	if (indexList.size() > 0)
		setDocumentation(m_model->data(indexList[0], Qt::ToolTipRole).toString());
}

void RenderSettingsDialog::update() {
	int index = ui->integratorBox->currentIndex();
	Properties integratorProps, samplerProps;
	bool needsUpdate = false;

	if (sender() == ui->samplerBox) {
		m_samplerNode->putProperties(samplerProps);
		needsUpdate = true;
	}

	if (sender() == ui->integratorBox) {
		m_integratorNode->putProperties(integratorProps);
		needsUpdate = true;
	}

	if (sender() == ui->rFilterBox ||
		sender() == ui->icBox ||
		sender() == ui->aiBox) {
		needsUpdate = true;
	}

	m_integratorNode = m_model->updateClass(m_integratorNode,
		ui->integratorBox->itemData(index).toList().at(0).toString(),
		ui->integratorBox->itemText(index));
	index = ui->samplerBox->currentIndex();
	m_samplerNode = m_model->updateClass(m_samplerNode,
		ui->samplerBox->itemData(index).toList().at(0).toString(),
		ui->samplerBox->itemText(index));
	index = ui->rFilterBox->currentIndex();
	m_rFilterNode = m_model->updateClass(m_rFilterNode,
		ui->rFilterBox->itemData(index).toList().at(0).toString(),
		ui->rFilterBox->itemText(index));

	if (ui->icBox->isChecked()) {
		m_icNode = m_model->updateClass(m_icNode,
			"IrradianceCacheIntegrator", tr("Irradiance Cache"));
	} else {
		m_icNode = m_model->updateClass(m_icNode, "", "");
	}

	if (ui->aiBox->isChecked()) {
		m_aiNode = m_model->updateClass(m_aiNode,
			"AdaptiveIntegrator", tr("Adaptive Integration"));
	} else {
		m_aiNode = m_model->updateClass(m_aiNode, "", "");
	}

	if (sender() == ui->integratorBox) {
		for (int i=0; i<m_integratorNode->childCount(); ++i) {
			TreeItem *treeItem = m_integratorNode->child(i);
			if (integratorProps.hasProperty(treeItem->getName().toStdString()))
				m_integratorNode->setProperty(treeItem->getName().toStdString(), integratorProps);
		}
	}

	if (sender() == ui->samplerBox) {
		for (int i=0; i<m_samplerNode->childCount(); ++i) {
			TreeItem *treeItem = m_samplerNode->child(i);
			if (samplerProps.hasProperty(treeItem->getName().toStdString()))
				m_samplerNode->setProperty(treeItem->getName().toStdString(), samplerProps);
		}
	}

	if (needsUpdate) {
		int row = 0;
		/* Make comboboxes etc editable by default */
		for (int i = 0; i < m_model->rowCount(); ++i) {
			QModelIndex index = m_model->index(i, 0);

			for (int j = 0; j < m_model->rowCount(index); ++j) {
				QModelIndex idx = m_model->index(j, 1, index);
				ui->treeView->openPersistentEditor(idx);
				QAbstractSpinBox *spinBox = qobject_cast<QAbstractSpinBox *>(ui->treeView->indexWidget(idx));
				if (spinBox) {
					QLineEdit *edit = spinBox->findChild<QLineEdit*>();
					if (row % 2 == 0)
						edit->setStyleSheet("background-color: palette(alternate-base);");
					edit->deselect();
				}
				row++;
			}
		}
	}
	ui->treeView->expandAll();

	dataChanged();
}

bool RenderSettingsDialog::resolutionHasChanged() const {
	return ui->resolutionBox->currentText() != m_originalResolution;
}

void RenderSettingsDialog::refresh() {
	bool valid = true;
	int pos;

	QString resolutionString(ui->resolutionBox->currentText());
	valid &= ui->resolutionBox->validator()->validate(resolutionString,pos)
		== QValidator::Acceptable;

	ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(valid);
}

void RenderSettingsDialog::load(const SceneContext *ctx) {
	const Scene *scene = ctx->scene.get();
	const Film *film = scene->getFilm();
	const Properties &rFilterProps = film->getReconstructionFilter()->getProperties();
	const Properties &samplerProps = scene->getSampler()->getProperties();
	const Integrator *integrator = scene->getIntegrator();
	Properties integratorProps = integrator->getProperties();

	if (integratorProps.getPluginName() == "adaptive") {
		ui->aiBox->setChecked(true);
		m_model->setProperties(m_aiNode, integratorProps);
		integrator = integrator->getSubIntegrator(0);
		integratorProps = integrator->getProperties();
	}

	if (integratorProps.getPluginName() == "irrcache") {
		ui->icBox->setChecked(true);
		m_model->setProperties(m_icNode, integratorProps);
		integrator = integrator->getSubIntegrator(0);
		integratorProps = integrator->getProperties();
	}

	ui->resolutionBox->lineEdit()->setText(QString("%1x%2")
		.arg(film->getCropSize().x).arg(film->getCropSize().y));
	m_originalResolution = ui->resolutionBox->lineEdit()->text();

	setComboBox(ui->integratorBox, integratorProps.getPluginName());
	setComboBox(ui->rFilterBox, rFilterProps.getPluginName());
	setComboBox(ui->samplerBox, samplerProps.getPluginName());
	update();

	m_model->setProperties(m_rFilterNode, rFilterProps);
	m_model->setProperties(m_samplerNode, samplerProps);
	m_model->setProperties(m_integratorNode, integratorProps);

	/* Make comboboxes etc editable by default */
	int row = 0;
	for (int i = 0; i < m_model->rowCount(); ++i) {
		QModelIndex index = m_model->index(i, 0);

		for (int j = 0; j < m_model->rowCount(index); ++j) {
			QModelIndex idx = m_model->index(j, 1, index);
			ui->treeView->openPersistentEditor(idx);
			QAbstractSpinBox *spinBox = qobject_cast<QAbstractSpinBox *>(ui->treeView->indexWidget(idx));
			if (spinBox) {
				QLineEdit *edit = spinBox->findChild<QLineEdit*>();
				if (row % 2 == 0)
					edit->setStyleSheet("background-color: palette(alternate-base);");
				edit->deselect();
			}
			row++;
		}
	}

	ui->treeView->expandAll();
}

void RenderSettingsDialog::apply(SceneContext *ctx) {
	Scene *scene = new Scene(ctx->scene);
	ref<Sensor> oldSensor = scene->getSensor();
	Film *oldFilm = oldSensor->getFilm();
	Properties filmProps = oldSensor->getFilm()->getProperties();
	ref<PluginManager> pluginMgr = PluginManager::getInstance();

	/* Temporarily set up a new file resolver */
	ref<Thread> thread = Thread::getThread();
	ref<FileResolver> oldResolver = thread->getFileResolver();
	ref<FileResolver> newResolver = oldResolver->clone();
	newResolver->prependPath(fs::absolute(scene->getSourceFile()).parent_path());
	thread->setFileResolver(newResolver);

	/* Configure the reconstruction filter */
	Properties rFilterProps(getPluginName(ui->rFilterBox));
	if (m_rFilterNode != NULL)
		m_rFilterNode->putProperties(rFilterProps);
	ref<ReconstructionFilter> rFilter = static_cast<ReconstructionFilter *>
		(pluginMgr->createObject(MTS_CLASS(ReconstructionFilter), rFilterProps));
	rFilter->configure();

	/* Configure the sampler */
	Properties samplerProps(getPluginName(ui->samplerBox));
	if (m_samplerNode != NULL)
		m_samplerNode->putProperties(samplerProps);
	ref<Sampler> sampler = static_cast<Sampler *>
		(pluginMgr->createObject(MTS_CLASS(Sampler), samplerProps));
	sampler->configure();

	/* Configure the integrator */
	Properties integratorProps(getPluginName(ui->integratorBox));
	if (m_integratorNode != NULL)
		m_integratorNode->putProperties(integratorProps);
	ref<Integrator> integrator = static_cast<Integrator *>
		(pluginMgr->createObject(MTS_CLASS(Integrator), integratorProps));
	integrator->configure();

	if (ui->icBox->isChecked()) {
		Properties icProps("irrcache");
		if (m_icNode != NULL)
			m_icNode->putProperties(icProps);
		ref<Integrator> ic = static_cast<Integrator *>
			(pluginMgr->createObject(MTS_CLASS(Integrator), icProps));
		ic->addChild(integrator);
		ic->configure();
		integrator = ic;
	}

	if (ui->aiBox->isChecked()) {
		Properties aiProps("adaptive");
		if (m_aiNode != NULL)
			m_aiNode->putProperties(aiProps);
		ref<Integrator> ai = static_cast<Integrator *>
			(pluginMgr->createObject(MTS_CLASS(Integrator), aiProps));
		ai->addChild(integrator);
		ai->configure();
		integrator = ai;
	}

	QStringList resolution = ui->resolutionBox->currentText().split('x');
	SAssert(resolution.size() == 2);
	Vector2i cropSize(
		std::max(1, resolution[0].toInt()),
		std::max(1, resolution[1].toInt()));

	/* Configure the film */
	Vector2i oldSize = oldFilm->getSize();
	Vector2i oldCropSize = oldFilm->getCropSize();
	Point2i oldCropOffset = oldFilm->getCropOffset();

	Vector2i size(math::roundToInt((oldSize.x * cropSize.x / (Float) oldCropSize.x)),
			      math::roundToInt((oldSize.y * cropSize.y / (Float) oldCropSize.y)));

	Point2i cropOffset(math::roundToInt((oldCropOffset.x * cropSize.x / (Float) oldCropSize.x)),
			           math::roundToInt((oldCropOffset.y * cropSize.y / (Float) oldCropSize.y)));

	filmProps.setInteger("width", size.x, false);
	filmProps.setInteger("height", size.y, false);

	if (size.x != cropSize.x || size.y != cropSize.y || cropOffset.x != 0 || cropOffset.y != 0) {
		filmProps.setInteger("cropWidth", cropSize.x, false);
		filmProps.setInteger("cropHeight", cropSize.y, false);
		filmProps.setInteger("cropOffsetX", cropOffset.x, false);
		filmProps.setInteger("cropOffsetY", cropOffset.y, false);
	} else {
		filmProps.removeProperty("cropWidth");
		filmProps.removeProperty("cropHeight");
		filmProps.removeProperty("cropOffsetX");
		filmProps.removeProperty("cropOffsetY");
	}

	ctx->originalSize = cropSize;

	ref<Film> film = static_cast<Film *> (pluginMgr->createObject(
			MTS_CLASS(Film), filmProps));
	film->addChild(rFilter);
	film->configure();

	if (cropSize.x != ctx->framebuffer->getWidth() ||
		cropSize.y != ctx->framebuffer->getHeight()) {
		ctx->framebuffer = new Bitmap(Bitmap::ERGBA, Bitmap::EFloat32, cropSize);
		ctx->framebuffer->clear();
		ctx->mode = EPreview;
	}

	/* Configure the sensor */
	Properties sensorProps = oldSensor->getProperties();

	if (oldSensor->getClass()->derivesFrom(MTS_CLASS(PerspectiveCamera))) {
		sensorProps.removeProperty("focalLength");
		sensorProps.setString("fovAxis", "y", false);
		sensorProps.setFloat("fov",
			static_cast<const PerspectiveCamera *>(oldSensor.get())->getYFov(), false);
	}

	ref<Sensor> newSensor = static_cast<Sensor *>
		(pluginMgr->createObject(MTS_CLASS(Sensor), sensorProps));
	newSensor->addChild(sampler);
	newSensor->addChild(film);
	newSensor->setMedium(oldSensor->getMedium());
	newSensor->configure();

	/* Update the scene with the newly constructed elements */
	scene->removeSensor(oldSensor);
	scene->addSensor(newSensor);
	scene->setSensor(newSensor);

	scene->setSampler(sampler);
	scene->setIntegrator(integrator);
	scene->configure();

	ctx->scene = scene;
	thread->setFileResolver(oldResolver);
}

RenderSettingsDialog::~RenderSettingsDialog() {
	delete ui;
}

void RenderSettingsDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
		ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

/* ====================== PropertyDelegate impl ====================== */

PropertyDelegate::PropertyDelegate(QObject *parent) : QStyledItemDelegate(parent) {
}

PropertyDelegate::~PropertyDelegate() {
}

QString	PropertyDelegate::displayText(const QVariant &value, const QLocale &locale) const {
	if (value.type() == QVariant::Bool) {
		#if defined(BOOLEAN_AS_COMBOBOXES)
			return value.toBool() ? tr("Yes") : tr("No");
		#else
			return QString("");
		#endif
	}
	return QStyledItemDelegate::displayText(value, locale);
}

void PropertyDelegate::updateWidgetData() {
	emit commitData((QWidget *) sender());
}

QWidget *PropertyDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option,
		const QModelIndex &index) const {
	if (index.data().type() == QVariant::Bool) {
		#if defined(BOOLEAN_AS_COMBOBOXES)
			QComboBox *cbox = new QComboBox(parent);
			cbox->addItem(tr("No"));
			cbox->addItem(tr("Yes"));
			return cbox;
		#else
			QCheckBox *box = new QCheckBox(parent);
			connect(box, SIGNAL(toggled(bool)), this, SLOT(updateWidgetData()));
			return box;
		#endif
	}

	QWidget *widget;
	if (index.data().type() == QVariant::Int)
		widget = new BetterSpinBox(parent);
	else if (index.data().type() == QVariant::Double)
		widget = new BetterDoubleSpinBox(parent);
	else
		widget = QStyledItemDelegate::createEditor(parent, option, index);

	#if defined(__OSX__)
		/* Don't draw focus halos on OSX, they're really distracting */
		if (widget != NULL && widget->testAttribute(Qt::WA_MacShowFocusRect))
			widget->setAttribute(Qt::WA_MacShowFocusRect, false);
		if (index.data().type() != QVariant::Bool) {
			widget->setAttribute(Qt::WA_MacMiniSize, true);
			widget->setStyleSheet("font-size: 13pt;");
		}
	#endif
	return widget;
}

void PropertyDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const {
	if (index.data().type() == QVariant::Bool) {
		#if defined(BOOLEAN_AS_COMBOBOXES)
			QComboBox *cbox = static_cast<QComboBox *>(editor);
			cbox->setCurrentIndex(index.data().toBool() ? 1 : 0);
		#else
			QCheckBox *cbox = static_cast<QCheckBox *>(editor);
			cbox->setChecked(index.data().toBool());
		#endif
		return;
	}

	QStyledItemDelegate::setEditorData(editor, index);
}

void PropertyDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
	const QModelIndex &index) const {
	if (index.data().type() == QVariant::Bool) {
		#if defined(BOOLEAN_AS_COMBOBOXES)
			QComboBox *cbox = static_cast<QComboBox *>(editor);
			model->setData(index, QVariant(cbox->currentIndex() == 1), Qt::EditRole);
		#else
			QCheckBox *cbox = static_cast<QCheckBox *>(editor);
			model->setData(index, QVariant(cbox->isChecked()), Qt::EditRole);
		#endif
		return;
	}
	QStyledItemDelegate::setModelData(editor, model, index);
}

void PropertyDelegate::updateEditorGeometry(QWidget *editor,
	const QStyleOptionViewItem &option, const QModelIndex &index) const {
	if (index.data().type() == QVariant::Bool) {
		editor->setGeometry(option.rect);
		return;
	}

	QStyledItemDelegate::updateEditorGeometry(editor, option, index);
}

QStringList RenderSettingsDialog::validateConfiguration() const {
	/* Ad-hoc verification until we have something better
	   (preferably specifiable by the plugins themselves) */
	QStringList messages;
	std::string integratorName = getPluginName(ui->integratorBox);
	std::string samplerName = getPluginName(ui->samplerBox);
	Properties integratorProps, samplerProps;
	m_integratorNode->putProperties(integratorProps);
	m_samplerNode->putProperties(samplerProps);

	if (samplerName != "independent") {
		if (integratorName == "pssmlt" || integratorName == "mlt")
			messages << "Error: Metropolis Light Transport-type algorithms only work with the independent sampler.";
		if (ui->aiBox->isChecked())
			messages << "Error: Adaptive integration requires the independent sampler.";
	}

	if ((samplerName == "ldsampler" || samplerName == "stratified") && integratorName == "ptracer")
		messages << "Error: the particle tracer does not support the stratified or low-discrepancy samplers!";

	if (samplerName == "halton" || samplerName == "hammersley") {
		if (integratorName == "bdpt")
			messages << "Error: the Bidirectional Path Tracer should not be used with the Halton/Hammersley samplers!";
		else if (integratorName == "erpt")
			messages << "Error: the Energy Redistribution Path Tracer should not be used with the Halton/Hammersley samplers!";
	}

	if (samplerName == "hammersley") {
		if (integratorName == "photonmapper")
			messages << "Error: the Hammersley sampler cannot be used with the photon mapper. Try the Halton sampler instead.";
	}

	if (ui->icBox->isChecked()) {
		if (integratorName != "direct" && integratorName != "path" && integratorName != "volpath" && integratorName != "volpath_simple" && integratorName != "photonmapper")
			messages << "Error: Irradiance Caching is not compatible with the selected integrator.";

	}
	if (ui->aiBox->isChecked()) {
		if (integratorName != "direct" && integratorName != "path" && integratorName != "volpath" && integratorName != "volpath_simple" && integratorName != "photonmapper")
			messages << "Error: Adaptive integration is not compatible with the selected integrator.";
	}
	if (integratorName == "ppm") {
		if ((samplerName == "independent" || samplerName == "ldsampler") && samplerProps.hasProperty("sampleCount")) {
			if (samplerProps.getInteger("sampleCount") > 4)
				messages << "Warning: are you sure you need more than 4 samples/pixel for progressive photon mapping? This will be slow..";
		} else if (samplerName == "stratified" && samplerProps.hasProperty("resolution")) {
			if (samplerProps.getInteger("resolution") > 2)
				messages << "Warning: are you sure you need more than 4 samples/pixel for progressive photon mapping? This will be slow..";
		}
	}
	return messages;
}

