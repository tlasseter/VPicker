package com.Sts.UI;

import com.Sts.DBTypes.StsClass;
import com.Sts.DBTypes.StsMainObject;
import com.Sts.DBTypes.StsMonitor;
import com.Sts.DBTypes.StsObject;
import com.Sts.MVC.Main;
import com.Sts.MVC.StsModel;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Table.StsDeleteRowNotifyListener;
import com.Sts.UI.Table.StsEditableTableModel;
import com.Sts.UI.Table.StsSelectRowNotifyListener;
import com.Sts.UI.Table.StsTableModelListener;
import com.Sts.Utilities.StsMath;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsObjectListPanel extends StsJPanel implements ActionListener, ListSelectionListener, StsSelectRowNotifyListener, StsDeleteRowNotifyListener, StsTableModelListener
{
    private StsModel model = null;
    private int[] selectedIndices = null;

    private ArrayList objectList = new ArrayList();
    StsGroupBox box = new StsGroupBox();
    JButton enableButton = new JButton();
    JButton propertiesBtn = new JButton();
    JButton exportBtn = new JButton();
    JButton disableButton = new JButton();
    StsTablePanelNew objectListTablePanel = null;

    StsMainObject[] objects;
    StsClass objectClass = null;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsObjectListPanel(StsModel model, StsClass objectClass)
    {
        try
        {
            this.model = model;
            this.objectClass = objectClass;
            buildTablePanel();
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void buildTablePanel()
    {
        String[] columnNames = {"isVisible", "name", "numberOfElements"};
        String[] columnTitles = {"Enabled", "Name", "Points"};
        StsEditableTableModel tableModel = new StsEditableTableModel(objectClass.getObjectClass(), columnNames, columnTitles, false);
        objectListTablePanel = new StsTablePanelNew(tableModel);
        objectListTablePanel.setSize(400, 100);
        objectListTablePanel.initialize();
        objectListTablePanel.setLabel("Objects");
    }

    public boolean initialize()
    {
        objectList.clear();
        objects = (StsMainObject[])objectClass.getCastObjectList();
        if(objects.length == 0)
            return false;
        int nObjects = objects.length;
        for(int n = 0; n < nObjects; n++)
            objectList.add(objects[n]);

        objectListTablePanel.replaceRows(objectList);
        objectListTablePanel.setSelectionIndex(0);
        objectListTablePanel.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        return true;
    }

    public StsMonitor getSelectedMonitor()
    {
        return null;
    }

    void jbInit() throws Exception
    {
        gbc.fill = gbc.BOTH;
        gbc.anchor = gbc.WEST;
        gbc.weighty = 1.0;
        gbc.fill = gbc.BOTH;
        gbc.gridwidth = 3;
        add(objectListTablePanel);

        objectListTablePanel.addSelectRowNotifyListener(this);
        objectListTablePanel.addTableModelListener(this);
        objectListTablePanel.addDeleteRowNotifyListener(this);

        enableButton.setText("Enable");
        propertiesBtn.setText("Properties");
        exportBtn.setText("Export");
        disableButton.setText("Disable");

        gbc.fill = gbc.NONE;
        gbc.anchor = gbc.CENTER;
        gbc.weighty = 0.0f;
        gbc.gridwidth = 1;
        box.addToRow(enableButton);
        if(Main.viewerOnly)
        {
            box.addToRow(propertiesBtn);
            propertiesBtn.addActionListener(this);
        }
        box.addToRow(exportBtn);
        box.addEndRow(disableButton);

        disableButton.addActionListener(this);
        exportBtn.addActionListener(this);
        enableButton.addActionListener(this);

        gbc.fill = gbc.HORIZONTAL;
        add(box);
    }

    public void actionPerformed(ActionEvent evt) {
        Object source = evt.getSource();
        if(source == disableButton)
        {
        	for(int i=0; i<selectedIndices.length; i++)
                objects[selectedIndices[i]].setIsVisible(false);
            initialize();
        }
        else if(source == propertiesBtn)
        {
            objects[selectedIndices[0]].popupPropertyPanel();
        }
        else if(source == enableButton)
        {
        	for(int i=0; i<selectedIndices.length; i++)
                objects[selectedIndices[i]].setIsVisible(true);
            initialize();
        }
        else if(source == exportBtn)
        {
        	for(int i=0; i<selectedIndices.length; i++)
                objects[selectedIndices[i]].export();
        }
    }

    public void valueChanged(ListSelectionEvent e)
	{
    	return;
 	}

    public void updateTable()
    {
    	if(objects == null) return;

    	int[] indices = new int[objects.length];
    	int cnt = 0;
    	int idx = 0;

        Iterator iter = objectList.iterator();
        while (iter.hasNext())
    	{
        	StsObject object = (StsMonitor) iter.next();

            objectListTablePanel.resetHighlight();
        	//monitorStatusTablePanel.setValueAt(monitor.getNumPolls(), monitor, "numPolls");

        	indices[cnt++] = idx;
        	idx++;
    	}
        StsMath.trimArray(indices,cnt);
    	repaint();
    }

    public void rowsSelected( int[] selected)
    {
    	selectedIndices = selected;
        return;
    }

    public void removeRows( int firstRow, int lastRow)
    {
        return;
    }

    public void deleteRow(int selectedRow)
    {
        return;
    }
}