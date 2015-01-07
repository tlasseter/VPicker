package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Table.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsWatchMonitorPanel extends StsJPanel implements ActionListener, ListSelectionListener,
StsSelectRowNotifyListener, StsDeleteRowNotifyListener, StsTableModelListener
{
    private JDialog parent = null;
    private StsModel model = null;
    private StsMonitor selectedMonitor = null;
    transient private boolean hasNewButton = false;

    private ArrayList monitorList = new ArrayList();
    JButton enableButton = new JButton();
    JButton editButton = new JButton();
    JButton newButton = new JButton();
    JButton disableButton = new JButton();
    StsTablePanelNew monitorStatusTablePanel = null;

    StsMonitor[] monitors;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsWatchMonitorPanel(StsModel model, JDialog parent, boolean newButton)
    {
        try
        {
            this.model = model;
            this.parent = parent;
            hasNewButton = newButton;
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
        String[] columnNames = {"name", "monitorObject", "source", "lastPollTimeString", "numPolls", "acceptedChanges", "rejectedChanges"};
        String[] columnTitles = {"Name", "Object", "Source", "Last Poll", "# Polls", "# Accepted", "# Rejected"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsMonitor.class, columnNames, columnTitles, false);
        monitorStatusTablePanel = new StsTablePanelNew(tableModel);
        monitorStatusTablePanel.setSize(400, 100);
        monitorStatusTablePanel.initialize();       
        monitorStatusTablePanel.setLabel("Monitor Objects");
    }

    public void initialize()
    {
        monitorList.clear();
        monitors = (StsMonitor[])model.getCastObjectList(StsMonitor.class);
        int nMonitors = monitors.length;
        for(int n = 0; n < nMonitors; n++)
            monitorList.add(monitors[n]);
        monitorStatusTablePanel.replaceRows(monitorList);
        monitorStatusTablePanel.setSelectionIndex(0);
        monitorStatusTablePanel.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
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
        gbc.gridwidth = 4;
        add(monitorStatusTablePanel, 1, 1);
        monitorStatusTablePanel.addSelectRowNotifyListener(this);
        monitorStatusTablePanel.addTableModelListener(this);
        monitorStatusTablePanel.addDeleteRowNotifyListener(this);

        enableButton.setText("Enable");
        disableButton.setText("Disable");
        editButton.setText("Edit...");
        newButton.setText("New...");
        gbc.fill = gbc.NONE;
        gbc.anchor = gbc.CENTER;
        gbc.weighty = 0.0f;
        gbc.gridwidth = 1;
        if(hasNewButton)
            addToRow(newButton);
        if(Main.viewerOnly)              // Otherwise use the Object Panel to edit monitors
            addToRow(editButton);
        addToRow(enableButton);
        addEndRow(disableButton);
        disableButton.addActionListener(this);
        enableButton.addActionListener(this);
        newButton.addActionListener(this);
        editButton.addActionListener(this);   
    }

    public void actionPerformed(ActionEvent evt) {
        Object source = evt.getSource();
        if(source == disableButton)
        	selectedMonitor.setEnable(false);
        else if(source == enableButton)
        	selectedMonitor.setEnable(true);
        else if(source == newButton)
        {
            model.win3d.getActionManager().launchWizard("com.Sts.Actions.Wizards.Monitor.StsMonitorWizard", "Monitors");
            parent.setVisible(false);
        }
        else if(source == editButton)
        {
            selectedMonitor.popupPropertyPanel();
            parent.setVisible(false);
        }
        selectedMonitor.updatePanel(selectedMonitor.lastPollTime);
    }

    public void valueChanged(ListSelectionEvent e)
	{
    	return;
 	}

    public void updateTable()
    {
    	if(monitors == null) return;

    	int[] indices = new int[monitors.length];
    	int cnt = 0;
    	int idx = 0;

        Iterator iter = monitorList.iterator();
        while (iter.hasNext())
    	{
        	StsMonitor monitor = (StsMonitor) iter.next();
        	if(monitor.hasChanged())
        	{

            	monitorStatusTablePanel.resetHighlight();
        		monitorStatusTablePanel.setValueAt(monitor.getLastPollTimeString(), monitor, "lastPollTimeString");
        		monitorStatusTablePanel.setValueAt(monitor.getNumPolls(), monitor, "numPolls");
        		monitorStatusTablePanel.setValueAt(monitor.getAcceptedChanges(), monitor, "acceptedChanges");
        		monitorStatusTablePanel.setValueAt(monitor.getRejectedChanges(), monitor, "rejectedChanges");
        		indices[cnt++] = idx;
        		monitor.resetChanged();
        	}
        	idx++;
    	}
        StsMath.trimArray(indices,cnt);
        //monitorStatusTablePanel.highlightRows(indices, Color.LIGHT_GRAY);
    	repaint();
    }

    public void rowsSelected( int[] selectedIndices)
    {
    	selectedMonitor = monitors[selectedIndices[0]];
        return;
    }

    public void removeRows( int firstRow, int lastRow)
    {
        return;
    }

    public void deleteRow(int selectedRow)
    {
    	selectedMonitor.delete();
        return;
    }
}
