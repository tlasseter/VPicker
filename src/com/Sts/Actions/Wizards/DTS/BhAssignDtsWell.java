package com.Sts.Actions.Wizards.DTS;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

/**
 * A list of DTS datasets that have been loaded are on the left.  The currently available wells are on the right.
 * By selecting a dataset on the left and a well on the right, the well will be associated with the dataset.
 */
public class BhAssignDtsWell extends StsWizardStep implements ActionListener, ListSelectionListener
{
    StsHeaderPanel header = new StsHeaderPanel();
    JPanel mainPanel = new JPanel(new BorderLayout());
    StsDoubleListPanel panel = new StsDoubleListPanel();

    public BhAssignDtsWell(StsWizard wizard)
    {
        super(wizard);
        setPanels(mainPanel, header);
        mainPanel.add(panel);
        panel.setLeftTitle("DTS Data");
        panel.setRightTitle("Well");
        panel.addActionListener(this);
        panel.addListSelectionListener(this);
        panel.getLeftButton().setToolTipText("Assign a well set to a dataset.");
        panel.getRightButton().setToolTipText("Remove a well association with a dataset.");

        header = (StsHeaderPanel)getHdrContainer();
        header.setTitle("DTS - Well Association");
        header.setSubtitle("Assign DTS Dataset to Well");
        header.setLink("http://www.s2ssystems.com/Marketing/AppLinks.html#DTS");        
        header.setInfoText(wizardDialog, "(1) Select a DTS data set from left list.\n" +
                "(2) Select the related well from the right list.\n" +
                "(3) Press the > button to assign the well to the DTS data set or\n" +
                "    Press the < button to unassign the well from the selected DTS data set.\n" +
                "(4) Press the Next>> Button to complete the operation.");
    }

    private void initWellList()
    {
        StsObject[] wells = model.getObjectList(StsWell.class);
        StsObject[] liveWells = model.getObjectList(StsLiveWell.class);

        StsListObjectItem[] items = new StsListObjectItem[wells.length + liveWells.length];
        for(int i = 0; i < wells.length; i++)
        {
            StsWell well = (StsWell)wells[i];
            String wellName = well.getName();
            items[i] = new StsListObjectItem(i, wellName, well);
        }
        for(int i = 0; i < liveWells.length; i++)
        {
            StsLiveWell well = (StsLiveWell)liveWells[i];
            String wellName = well.getName();
            items[i] = new StsListObjectItem(i, wellName, well);
        }
        panel.setRightItems(items);
    }

    private void initDataSetList()
    {
        StsObject[] dataSets = ((BhDtsLoadWizard)wizard).getDatasets();
        StsListObjectItem[] items = new StsListObjectItem[dataSets.length];
        for(int i = 0; i < dataSets.length; i++)
        {
            StsTimeLogCurve dataset = (StsTimeLogCurve)dataSets[i];
            String dataSetName = dataset.getName();
            StsWell well = dataset.getWell();
            if(dataset.getWell() != null)
                dataSetName = dataSetName + "(" + well.getName() + ")";
            else
                dataSetName = dataSetName + "()";

            items[i] = new StsListObjectItem(i, dataSetName, dataset);
        }
        panel.setLeftItems(items);
    }

    public String trimName(String fullString)
    {
        int index = fullString.indexOf("(");
        if(index > 0)
            return fullString.substring(0, index);
        return null;
    }

    public boolean start()
    {
        refreshLists();
        enableFinish();
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public void refreshLists()
    {
        initDataSetList();
        initWellList();
        refreshButtons();
    }

    public void valueChanged(ListSelectionEvent e)
    {
        refreshButtons();
    }

    /** datasets are on the left and wells on the right. */
    public void refreshButtons()
    {
        Object[] datasetItems = panel.getSelectedLeftItems();
        Object[] wellItems = panel.getSelectedRightItems();
        int nDatasetItems = datasetItems == null ? 0 : datasetItems.length;
        int nWellItems = wellItems == null ? 0 : wellItems.length;
        panel.enableLeftButton(nDatasetItems > 0 && nDatasetItems == nWellItems);
        StsTimeLogCurve dataset = null;

        if(nDatasetItems > 0 && nWellItems == 0)
        {
            for(int i = 0; i < nDatasetItems; i++)
            {
                StsListObjectItem item = (StsListObjectItem)datasetItems[i];
                try
                {
                    dataset = (StsTimeLogCurve)item.object;
                }
                catch(Exception ex){ }
                if(dataset != null && dataset.getWell() != null)
                {
                    panel.enableRightButton(true);
                    return;
                }
            }
        }
        panel.enableRightButton(false);
    }

    public void actionPerformed(ActionEvent e)
    {
        if(e.getActionCommand().equals("<"))
        {
            Object[] datasetItems = panel.getSelectedLeftItems();
            Object[] wellItems = panel.getSelectedRightItems();
            int nDatasetItems = datasetItems == null ? 0 : datasetItems.length;
            int nWellItems = wellItems == null ? 0 : wellItems.length;
            if(nWellItems == nDatasetItems)
            {
                for(int i = 0; i < nDatasetItems; i++)
                {
                    StsListObjectItem item = (StsListObjectItem)datasetItems[i];
                    StsTimeLogCurve dataset = null;
                    try { dataset = (StsTimeLogCurve)item.object; }
                    catch(Exception ex){ }

                    item = (StsListObjectItem)wellItems[i];
                    StsWell well = (StsWell)item.object;

                    if(dataset != null && well != null)
                    {
                        // well.addLogCurve(dataset);
                        dataset.setWell(well);
                    }
                }
                refreshLists();
            }
        }
        else if(e.getActionCommand().equals(">"))
        {
            Object[] datasetItems = panel.getSelectedLeftItems();
            int nDatasetItems = datasetItems == null ? 0 : datasetItems.length;
            for(int i = 0; i < nDatasetItems; i++)
            {
                StsListObjectItem item = (StsListObjectItem)datasetItems[i];
                StsTimeLogCurve dataset = null;
                try
                { dataset = (StsTimeLogCurve)item.object; }
                catch(Exception ex){ }
                if(dataset != null)
                {
                    StsWell well = dataset.getWell();
                    if(well != null)
                    {
                        // dataset.setWell(null);
                        well.deleteLogCurve(dataset);
                    }
                }
            }
            refreshLists();
        }
    }
}
