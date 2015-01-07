package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectSensorPanel extends StsJPanel implements ActionListener
{
    private StsMonitorWizard wizard;
    private StsSelectSensor wizardStep;

    private StsModel model = null;
    private StsSensor selectedSensor = null;

    JList sensorList = new JList();
    DefaultListModel sensorListModel = new DefaultListModel();
    protected JScrollPane jScrollPane1;
    JButton newButton = new JButton();

    Object[] sensors;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsSelectSensorPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMonitorWizard)wizard;
        this.wizardStep = (StsSelectSensor)wizardStep;
        try
        {
            jbInit();
            //initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        sensorListModel.removeAllElements();

        sensors = (Object[])model.getObjectList(StsDynamicSensor.class);
        sensors = (Object[])StsMath.arrayAddArray(sensors, model.getObjectList(StsStaticSensor.class));

        int nElements = sensors.length;
        for(int n = 0; n < nElements; n++)
            sensorListModel.addElement(((StsSensor)sensors[n]).getName());

        sensorList.setModel(sensorListModel);
    }

    public StsSensor getSelectedSensor()
    {
        if(sensorList.isSelectionEmpty()) return null;
        return (StsSensor)sensors[sensorList.getSelectedIndex()];
    }

    void jbInit() throws Exception
    {
        this.gbc.fill = gbc.BOTH;
        this.gbc.anchor = gbc.WEST;
        this.gbc.weighty = 1.0;

        sensorList.setBorder(BorderFactory.createEtchedBorder());
        newButton.setText("New Sensor...");
        jScrollPane1 = new JScrollPane(sensorList);
        jScrollPane1.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        jScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        addEndRow(jScrollPane1);

        this.gbc.fill = gbc.NONE;
        this.gbc.anchor = gbc.CENTER;
        this.gbc.weighty = 0.0f;
        addEndRow(newButton);
        newButton.addActionListener(this);

        //newButton.setEnabled(false);   // Not fully tested yet."
    }

    public void actionPerformed(ActionEvent evt) {
        Object source = evt.getSource();
        if (source == newButton)
        	wizard.newSensor();
    }
}