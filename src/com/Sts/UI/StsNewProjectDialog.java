package com.Sts.UI;

import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsNewProjectDialog extends JDialog
{
    StsJPanel panel = StsJPanel.addInsets();
    StsButton okBtn = new StsButton("  OK  ", "Accept these units.", this, "executeOK");
    StsButton cancelBtn = new StsButton(" Cancel ", "Accept default units.", this, "executeCancel");
    StsGroupBox unitsGroupBox = new StsGroupBox("Project Units");
    StsGroupBox colorsGroupBox = new StsGroupBox("Default Colors");
	StsGroupBox buttonsGroupBox = new StsGroupBox();

    byte timeUnits, depthUnits, xyUnits;
    float depthDatum = 0.0f, timeDatum = 0.0f;
//	Color defaultWellColor = Color.RED;

    StsProject project = null;

    StsComboBoxFieldBean xyUnitsBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean timeUnitsBean = new StsComboBoxFieldBean();
    StsComboBoxFieldBean depthUnitsBean = new StsComboBoxFieldBean();

    StsFloatFieldBean depthDatumBean = new StsFloatFieldBean();
    StsFloatFieldBean timeDatumBean = new StsFloatFieldBean();

    //StsComboBoxFieldBean initialDomainBean = new StsComboBoxFieldBean();

	StsColorComboBoxFieldBean wellColorComboBoxBean = new StsColorComboBoxFieldBean();

    boolean canceled = false;

    public StsNewProjectDialog(StsProject project)
    {
        super((Frame)null, "New Project Setup", true);
        try
        {
            this.project = project;
            timeUnits = project.getTimeUnits();
            xyUnits = project.getXyUnits();
            depthUnits = project.getDepthUnits();

            //initialDomainBean.classInitialize(project, "zDomainString", "Initial Domain", StsParameters.TD_STRINGS);
            timeUnitsBean.initialize(this, "timeUnits", " Time Units:  ", StsParameters.TIME_STRINGS);
            xyUnitsBean.initialize(this, "xyUnits", "   XY Units:  ", StsParameters.DIST_STRINGS);
            depthUnitsBean.initialize(this, "depthUnits", "Depth Units:  ", StsParameters.DIST_STRINGS);
	//		wellColorListBean.initializeColors(this, "defaultWellColor", "Default Well Color", colors);

            depthDatumBean.initialize(this,"depthDatum",true, "Depth Datum:  ",true);
            depthDatumBean.setRangeFixStep(-1000.0f,1000.0f,10.0f);
            timeDatumBean.initialize(this,"timeDatum",true, "Time Datum:  ",true);
            timeDatumBean.setRangeFixStep(-1000.0f,1000.0f,10.0f);

            jbInit();
			setSize(200, 400);
            StsToolkit.centerComponentOnScreen(this);
            pack();

            timeUnitsBean.setSelectedItem(project.getTimeUnitString());
            xyUnitsBean.setSelectedItem(project.getXyUnitString());
            depthUnitsBean.setSelectedItem(project.getDepthUnitString());
            //initialDomainBean.setSelectedItem(project.getStringForZDomain());
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        //unitsGroupBox.add(initialDomainBean);
        unitsGroupBox.add(depthUnitsBean);
        //unitsGroupBox.add(depthDatumBean);
		unitsGroupBox.add(xyUnitsBean);
		unitsGroupBox.add(timeUnitsBean);
        //unitsGroupBox.add(timeDatumBean);

//	    colorsGroupBox.add(wellColorListBean);

//	    buttonsGroupBox.gbc.fill = GridBagConstraints.HORIZONTAL;
	    buttonsGroupBox.addToRow(okBtn);
		buttonsGroupBox.addEndRow(cancelBtn);

		panel.gbc.fill = GridBagConstraints.HORIZONTAL;
	    panel.add(unitsGroupBox);
//		panel.add(colorsGroupBox);
		panel.add(buttonsGroupBox);

        getContentPane().add(panel);
    }

	public void executeOK()
	{
		project.setTimeUnits(timeUnits);
		project.setXyUnits(xyUnits);
		project.setDepthUnits(depthUnits);
        //project.setTimeDatum(timeDatum);
        //project.setDepthDatum(depthDatum);
//		project.setDefaultWellColor(defaultWellColor);
		setVisible(false);
	}

	public void executeCancel()
	{
		canceled = true;
		setVisible(false);
	}
    public void setTimeDatum(float datum) { timeDatum = datum; }
    public float getTimeDatum() { return timeDatum; }
    public void setDepthDatum(float datum) { depthDatum = datum; }
    public float getDepthDatum() { return depthDatum; }
    public void setTimeUnits(String units) { timeUnits = StsParameters.getTimeUnitsFromString(units); }
    public String getTimeUnits() { return StsParameters.getTimeUnitString(timeUnits); }
    public void setXyUnits(String units) { xyUnits = StsParameters.getDistanceUnitsFromString(units); }
    public String getXyUnits() { return StsParameters.getDistanceUnitString(xyUnits); }
    public void setDepthUnits(String units) { depthUnits = StsParameters.getDistanceUnitsFromString(units); }
    public String getDepthUnits() { return StsParameters.getDistanceUnitString(depthUnits); }
//	public void setDefaultWellColor(Color color) { this.defaultWellColor = color; }
//	public Color getDefaultWellColor() { return defaultWellColor; }

    public boolean wasCanceled() { return canceled; }

	static public void main(String[] args)
	{
		StsProject project = new StsProject();
		StsNewProjectDialog dialog = new StsNewProjectDialog(project);
		dialog.pack();
		dialog.setVisible(true);
	}
}
