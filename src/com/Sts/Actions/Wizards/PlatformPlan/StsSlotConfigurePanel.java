package com.Sts.Actions.Wizards.PlatformPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSlotConfigurePanel extends JPanel implements TableModelListener
{
    public StsPlatformPlanWizard wizard;
    private StsSlotConfigure wizardStep;

    private StsModel model = null;
    private StsPlatform platform = null;

    StsGroupBox slotConfigBox = new StsGroupBox("Slot Configuration");
    DefaultTableModel tableModel = new DefaultTableModel();
    private JTable slotTable = new JTable();
    private JScrollPane tablePane = new JScrollPane();

    String[][] cells = null;
    String[] colLabels = null;

    StsGroupBox slotGraphicBox = new StsGroupBox("Slot Graphic");
    StsPlatformConfigPanel slotGraphicPanel = new StsPlatformConfigPanel();

    public StsSlotConfigurePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsPlatformPlanWizard)wizard;
        this.wizardStep = (StsSlotConfigure)wizardStep;
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
        initializeTable();
    }

    void jbInit() throws Exception
    {
        // Configure some of JTable's paramters
        slotTable.setShowHorizontalLines(true);
        slotTable.setShowVerticalLines(true);
        slotTable.setCellSelectionEnabled(true);

        // Set selection parameters
        slotTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

        // Add the table to a scrolling pane
        tablePane = slotTable.createScrollPaneForTable( slotTable );
        tablePane.setPreferredSize(new Dimension(300, 200));

        slotConfigBox.add(tablePane, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
              , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 5, 5), 0, 0));


        slotGraphicPanel.setBackground(Color.WHITE);
        slotGraphicPanel.setSize(300,300);
        slotGraphicBox.add(slotGraphicPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
              , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 5, 5, 5), 300, 300));
        add(slotGraphicBox);
        add(slotConfigBox);
    }

    public void tableChanged(TableModelEvent e)
    {
        int rowIndex = slotTable.getSelectedRow();
        int colIndex = slotTable.getSelectedColumn();
        if((rowIndex > -1) && (colIndex > -1))
        {
            platform.setSlotName((rowIndex*platform.getNCols()) + colIndex,
                                 (String) slotTable.getModel().getValueAt(rowIndex, colIndex));
            platform.drawConfiguration();
        }
    }

    private void initializeTable()
    {
        int i, j;
        String[] slotNames = platform.getSlotNames();
        Object[] row = new Object[platform.getNCols()];

        tableModel = new DefaultTableModel();
        slotTable.setModel(tableModel);

        for (i = 0; i < platform.getNCols(); i++)
            tableModel.addColumn(String.valueOf(i + 1));

        tableModel.setNumRows(platform.getNRows());
        for(j=0; j< platform.getNRows(); j++)
        {
            for (i = 0; i < platform.getNCols(); i++)
            {
                slotTable.getModel().setValueAt(new String(slotNames[(j*platform.getNCols()) + i]), j, i);
            }
        }
        tableModel.addTableModelListener(this);
    }

}
