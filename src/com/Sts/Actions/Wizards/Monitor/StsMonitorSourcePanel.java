package com.Sts.Actions.Wizards.Monitor;

import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;

import javax.swing.*;
import java.awt.event.ActionListener;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.DBTypes.StsMonitor;
import com.Sts.DBTypes.StsSensor;
import com.Sts.DBTypes.StsDynamicSensor;
import com.Sts.DBTypes.StsStaticSensor;
import com.Sts.MVC.StsModel;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.Utilities.StsMath;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsMonitorSourcePanel extends StsJPanel implements ActionListener
{
    private StsMonitorWizard wizard;
    private StsMonitorSource wizardStep;

    private StsModel model = null;

    ButtonGroup typeGrp = new ButtonGroup();
    StsGroupBox typeBox = new StsGroupBox("Select Type");
    JRadioButton fromFileBtn = new JRadioButton("Define from File");
    JRadioButton manualBtn = new JRadioButton("Define Manually");

    ButtonGroup objectGrp = new ButtonGroup();
    StsGroupBox objectBox = new StsGroupBox("Select Data Type");
    JRadioButton sensorBtn = new JRadioButton("Sensor");
    JRadioButton wellLogBtn = new JRadioButton("Well & Logs");
    GridBagLayout gridBagLayout = new GridBagLayout();

    private boolean isManual = true;
    private byte monitorType = StsMonitor.SENSOR;

    public StsMonitorSourcePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMonitorWizard)wizard;
        this.wizardStep = (StsMonitorSource)wizardStep;
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
    }

    void jbInit() throws Exception
    {
        setLayout(gridBagLayout);

        typeBox.gbc.anchor = gbc.WEST;
        typeBox.gbc.fill = gbc.NONE;
        objectBox.gbc.anchor = gbc.WEST;
        objectBox.gbc.fill = gbc.NONE;

        typeBox.add(fromFileBtn);
        typeBox.add(manualBtn);
        objectBox.add(sensorBtn);
        objectBox.add(wellLogBtn);

        typeGrp.add(fromFileBtn);
        typeGrp.add(manualBtn);
        objectGrp.add(sensorBtn);
        objectGrp.add(wellLogBtn);

        fromFileBtn.addActionListener(this);
        manualBtn.addActionListener(this);
        sensorBtn.addActionListener(this);
        wellLogBtn.addActionListener(this);

        gbc.fill = gbc.HORIZONTAL;
        gbc.anchor = gbc.NORTH;
        addEndRow(typeBox);
        addEndRow(objectBox);

        fromFileBtn.setSelected(false);
        manualBtn.setSelected(true);

        sensorBtn.setSelected(true);
        wellLogBtn.setSelected(false);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == fromFileBtn)
            isManual = false;
        else if(source == manualBtn)
            isManual = true;
        else if(source == sensorBtn)
            monitorType = StsMonitor.SENSOR;
        else if(source == wellLogBtn)
            monitorType = StsMonitor.WELL;
    }

    public boolean isManual()
    {
       return isManual;
    }

    public byte getMonitorType()
    {
        return monitorType;
    }
}