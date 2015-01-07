package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.DBTypes.*;
import com.Sts.MVC.StsModel;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.Utilities.StsMath;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectWellPanel extends StsJPanel implements ActionListener
{
    private StsMonitorWizard wizard;
    private StsSelectWell wizardStep;

    private StsModel model = null;
    private StsSensor selectedWell = null;

    JList wellList = new JList();
    DefaultListModel wellListModel = new DefaultListModel();
    protected JScrollPane jScrollPane1;
    JButton newButton = new JButton();

    Object[] wells;

    public StsSelectWellPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsMonitorWizard)wizard;
        this.wizardStep = (StsSelectWell)wizardStep;
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
        wellListModel.removeAllElements();

        wells = (Object[])model.getObjectList(StsWell.class);
        //wells = (Object[])StsMath.arrayAddArray(wells, model.getObjectList(StsDrillingWells.class));

        int nElements = wells.length;
        for(int n = 0; n < nElements; n++)
            wellListModel.addElement(((StsWell)wells[n]).getName());

        wellList.setModel(wellListModel);
    }

    public StsWell getSelectedWell()
    {
        if(wellList.isSelectionEmpty()) return null;
        return (StsWell)wells[wellList.getSelectedIndex()];
    }

    void jbInit() throws Exception
    {
        this.gbc.fill = gbc.BOTH;
        this.gbc.anchor = gbc.WEST;
        this.gbc.weighty = 1.0;

        wellList.setBorder(BorderFactory.createEtchedBorder());
        newButton.setText("New Drilling Well...");
        jScrollPane1 = new JScrollPane(wellList);
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
        	wizard.newWell();
    }
}