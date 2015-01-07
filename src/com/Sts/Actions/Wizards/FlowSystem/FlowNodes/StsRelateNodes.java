
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.FlowSystem.FlowNodes;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

public class StsRelateNodes extends StsWizardStep implements ActionListener, ListSelectionListener
{
    StsHeaderPanel header = new StsHeaderPanel();
    StsSensor[] modelSurfaces = null;
	JPanel mainPanel = new JPanel(new BorderLayout());
	StsDoubleListPanel panel = new StsDoubleListPanel();

    public StsRelateNodes(StsWizard wizard)
    {
		super(wizard);
		setPanels(mainPanel, header);
        mainPanel.add(panel);
        panel.setLeftTitle("Tanks");
        panel.setRightTitle("Pumps");
        panel.addActionListener(this);
        panel.addListSelectionListener(this);
        panel.getLeftButton().setToolTipText("Create node between the selected pump & tank.");
        panel.getRightButton().setToolTipText("Remove node between the selected tank & pump.");

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Node Construction");
        header.setSubtitle("Assign Pumps to Tanks");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#FlowNodes");        
        header.setInfoText(wizardDialog,"(1) Select a tank from left list.\n" +
                                  "(2) Select the related pump from the right list.\n" +
                                  "(3) Press the > button to assign the pump to the tank or\n" +
                                  "    Press the < button to unassign the pump from the selected marker.\n" +
                                  "(4) Press the Next>> Button to complete node construction.");
    }

    private void initPumpList()
    {
        StsObject[] pumps = model.getObjectList(StsPump.class);
        StsListItem[] items = new StsListItem[pumps.length];
        for(int i=0; i<pumps.length; i++)
        {
        	StsPump pump = (StsPump)pumps[i];
        	String pumpName = pump.getName();
        	items[i] = new StsListItem(i, pumpName);
        }
        panel.setRightItems(items);
    }
    
    private void initTankList()
    {
        StsObject[] tanks = model.getObjectList(StsTank.class);
        StsListItem[] items = new StsListItem[tanks.length];
        for(int i=0; i<tanks.length; i++)
        {
        	StsTank tank = (StsTank)tanks[i];
        	String tankName = tanks[i].getName();
        	if(tank.getPump() != null)
        		tankName = tankName + "(" + tank.getPump().getName() + ")";
        	else
        		tankName = tankName + "()";
        	
        	items[i] = new StsListItem(i, tankName);
        }
        panel.setLeftItems(items);
    }
    
    public String trimName(String fullString)
    {
    	int index = fullString.indexOf("(");
        if( index > 0 )
			return fullString.substring(0, index);
        return null;
    }
    
    public boolean start()
    {
		refreshLists();
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public void refreshLists()
    {
    	initTankList();
    	initPumpList();
        refreshButtons();
    }

    public StsPump[] getSelectedPumps()
    {
        StsPump[] pumps = null;

        try
        {
            Object[] items = panel.getSelectedRightItems();
            int nItems = items == null ? 0 : items.length;
            if( nItems == 0) return null;

            pumps = new StsPump[nItems];
            for( int i=0; i<nItems; i++ )
            {
                StsListItem item = (StsListItem)items[i];
                pumps[i] = (StsPump)model.getObjectWithName(StsPump.class, item.name);
            }
    	    return pumps;
        }
        catch(Exception e)
        {
            StsException.outputException("StsRelateNodes.getSelectedPumps() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public void valueChanged(ListSelectionEvent e)
    {
    	refreshButtons();
    }

    public void refreshButtons()
    {
		Object[] tankItems = panel.getSelectedLeftItems();
		Object[] pumpItems = panel.getSelectedRightItems();
    	int nTankItems = tankItems == null ? 0 : tankItems.length;
    	int nPumpItems = pumpItems == null ? 0 : pumpItems.length;
		panel.enableLeftButton(nTankItems > 0 && nTankItems == nPumpItems);
		StsTank tank = null;
		
		if(nTankItems > 0 && nPumpItems == 0)
        {
			for( int i=0; i<nTankItems; i++ )
            {
            	StsListItem item = (StsListItem) tankItems[i];
                try 
                { 
                	tank = (StsTank)model.getObjectWithName(StsTank.class, trimName(item.name));
                }
                catch(Exception ex) { }
                if( tank != null && tank.getPump() != null )
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
    	if( e.getActionCommand().equals("<") )
        {
			Object[] tankItems = panel.getSelectedLeftItems();
			Object[] pumpItems = panel.getSelectedRightItems();
	    	int nTankItems = tankItems == null ? 0 : tankItems.length;
    		int nPumpItems = pumpItems == null ? 0 : pumpItems.length;
            if( nPumpItems == nTankItems )
            {
                for( int i=0; i<nTankItems; i++ )
                {
                    StsListItem item = (StsListItem) tankItems[i];
                    StsTank tank = null;
                    try { tank = (StsTank)model.getObjectWithName(StsTank.class, trimName(item.name)); }
                    catch(Exception ex) { }

                    item = (StsListItem) pumpItems[i];
                    StsPump pump = (StsPump)model.getObjectWithName(StsPump.class, item.name);

                    if( tank != null && pump != null )
                    {
                    	tank.setPump(pump);
                    }
                }
                refreshLists();
            }
        }
        else if( e.getActionCommand().equals(">") )
        {
			Object[] tankItems = panel.getSelectedLeftItems();
	    	int nTankItems = tankItems == null ? 0 : tankItems.length;
            for( int i=0; i<nTankItems; i++ )
            {
                StsListItem item = (StsListItem) tankItems[i];
                StsTank tank = null;
                try { tank = (StsTank)model.getObjectWithName(StsTank.class, trimName(item.name)); }
                catch(Exception ex) { }
                if( tank != null )
                {
                	StsPump oldPump = tank.getPump();
                    tank.setPump(null);
                }
            }
            refreshLists();
        }
    }
}
