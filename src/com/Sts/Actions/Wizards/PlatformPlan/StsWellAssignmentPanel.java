package com.Sts.Actions.Wizards.PlatformPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
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

public class StsWellAssignmentPanel extends JPanel implements ListSelectionListener, ActionListener
{
    public StsPlatformPlanWizard wizard;
    private StsWellAssignment wizardStep;

    private StsModel model = null;
    private StsPlatform platform = null;
    String[] wellnames = null;

    StsGroupBox slotConfigBox = new StsGroupBox("Well Assignment");
    private StsWell selectedWell = null;
    private int selectedIndex;
    private DefaultListModel wellsListModel = new DefaultListModel();
    JList wellsList = new JList(wellsListModel);
    private JButton clearAllBtn = new JButton("Clear All");
    private JButton clearCurrentSlotBtn = new JButton("Clear Current");

    private JScrollPane tablePane = new JScrollPane();

    StsGroupBox slotGraphicBox = new StsGroupBox("Slot Graphic");
    StsPlatformConfigPanel slotGraphicPanel = new StsPlatformConfigPanel();

    public StsWellAssignmentPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsPlatformPlanWizard)wizard;
        this.wizardStep = (StsWellAssignment)wizardStep;
        this.model = wizard.getModel();
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
        this.platform = wizard.getPlatform();
        platform.configurePlatform();

        slotGraphicPanel.setPlatform(platform);
        slotGraphicPanel.setLabelSlots(true);

        StsProject project = wizard.getModel().getProject();
        platform.setCanvas(slotGraphicPanel);
        platform.drawConfiguration();

        model = wizard.getModel();
        StsWell[] wells = (StsWell[])model.getCastObjectList(StsWell.class);
        wellnames = new String[wells.length];
        StsPlatformClass pc = (StsPlatformClass)model.getStsClass(StsPlatform.class);
        int count = 0;
        for(int i=0; i<wells.length; i++)
        {
            String wellname = wells[i].getName();
            StsPlatform wellPlatform = pc.getWellPlatform(wellname);
            if(wellPlatform == null || wellPlatform == platform)
                wellnames[count++] = wellname;
        }
        wellnames = (String[])StsMath.trimArray(wellnames, count);
        clearSlotAssignment();
    }

    void jbInit() throws Exception
    {
        tablePane.setPreferredSize(new Dimension(400, 200));
        wellsList.addListSelectionListener(this);

        tablePane.getViewport().add(wellsList, null);
        slotConfigBox.add(tablePane, new GridBagConstraints(0, 0, 2, 1, 1.0, 1.0
              , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 5, 5), 0, 0));
        slotConfigBox.add(clearAllBtn, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
              , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 5, 5, 5), 0, 0));
        slotConfigBox.add(clearCurrentSlotBtn, new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0
              , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 5, 5, 5), 0, 0));
        clearAllBtn.addActionListener(this);
        clearCurrentSlotBtn.addActionListener(this);

        slotGraphicPanel.setBackground(Color.WHITE);
        slotGraphicPanel.setSize(300,300);
        slotGraphicBox.add(slotGraphicPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
              , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 5, 5, 5), 300, 300));
        add(slotGraphicBox);
        add(slotConfigBox);
    }

    public void valueChanged(ListSelectionEvent e)
    {
        int i;
        Object source = e.getSource();
        String[] slotNames = platform.getSlotNames();

        wellsList.removeListSelectionListener(this);
        selectedIndex = wellsList.getSelectedIndex();
        if(platform.getCurrentSlotIndex() == -1)
        {
            new StsMessage(wizard.frame, StsMessage.WARNING, "Must select slot first.");
            wellsList.clearSelection();
            wellsList.addListSelectionListener(this);
            return;
        }

        if(platform.addWellWithSlotName(slotNames[platform.getCurrentSlotIndex()], wellnames[selectedIndex]))
            wellsListModel.setElementAt(wellnames[selectedIndex] + " in slot " + platform.getSlotName(platform.getCurrentSlotIndex()), selectedIndex);
        wellsList.addListSelectionListener(this);
        platform.drawConfiguration();
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        wellsList.removeListSelectionListener(this);
        if(source == clearAllBtn)
        {
            platform.clearSlotAssignment();
            clearSlotAssignment();
        }
        else if(source == clearCurrentSlotBtn)
        {
            wellsListModel.setElementAt(wellnames[selectedIndex], selectedIndex);
            platform.clearSlotAssignment(wellnames[selectedIndex]);
        }
        wellsList.addListSelectionListener(this);
    }

    public void clearSlotAssignment()
    {
        int nWells = wellnames.length;
        wellsListModel.removeAllElements();
        for(int i=0; i<nWells; i++)
        {
            if(platform.isWellOnPlatform(wellnames[i]))
                wellsListModel.addElement(wellnames[i] + " in slot " + platform.getSlotName(wellnames[i]));
            else
                wellsListModel.addElement(wellnames[i]);
        }
    }
}
