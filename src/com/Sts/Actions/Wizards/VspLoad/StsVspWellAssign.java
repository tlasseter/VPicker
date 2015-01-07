
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VspLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

public class StsVspWellAssign extends StsWizardStep	implements ActionListener, ListSelectionListener
{
    StsHeaderPanel header = new StsHeaderPanel();
	JPanel mainPanel = new JPanel(new BorderLayout());
    StsVsp[] vsps = null;
	StsDoubleListPanel panel = new StsDoubleListPanel();
    StsVspLoadWizard wizard = null;

    public StsVspWellAssign(StsVspLoadWizard wizard)
    {
        super(wizard);
        this.wizard = wizard;
        setPanels(mainPanel, header);
        mainPanel.add(panel);
        panel.setLeftListRenderer(new StsColorListRenderer());
        panel.setRightListRenderer(new StsColorListRenderer());
        panel.setLeftTitle("Selected VSP:");
        panel.setRightTitle("Available Wells:");
        panel.addActionListener(this);
        panel.addListSelectionListener(this);
        panel.getLeftButton().setToolTipText("Assign well to selected VSP");
        panel.getRightButton().setToolTipText("Remove well assignment from selected VSP");
        panel.setPreferredSize(new Dimension(400,300));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("VSP Load");
        header.setSubtitle("Assign Wells to VSP Data");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VspLoad");                        
        header.setInfoText(wizardDialog,"(1) Select a VSP from left list.\n" +
                                  "(2) Select the related well from the right list.\n" +
                                  "   ***** Well list is generated all wells loaded in project *****\n" +
                                  "   ***** Psuedo wells can be created prior to assignment *****\n" +
                                  "(3) Press create button to create vertical psuedo well\n" +
                                  "(4) Press the > button to assign the well to the selected VSP or\n" +
                                  "    Press the < button to unassign the well from the selected VSP.\n" +
                                  "(5) Press the Next>> Button to complete the VSP load.");
    }

    public void initVspList()
    {
        vsps = (StsVsp[]) model.getCreateStsClass(StsVsp.class).getChildren();
        if(vsps != null )
        {
            int nVsps = vsps == null ? 0 : vsps.length;
            if( nVsps > 0 )
            {
                StsColorListItem[] items = new StsColorListItem[nVsps];
                for( int i=0; i<nVsps; i++ )
                {
                    StsVsp vsp = (StsVsp)vsps[i];
                    String vspName = vsp.getName();
                    StsWell well = vsp.getWell();
                    if(well == null)
                    {
                        well = findWell(vsp, (int) 10);
                        if (well != null)
                        {
                            setVspDisplayFlags(well);
                            vsp.setWell(well);
                        }
                    }
                    if (well == null)
                    {
                        vspName += " ( )";
                        items[i] = new StsColorListItem(StsColor.RED, vspName);
                    }
                    else
                    {
                        vspName += " (" + well.getName() + ")";
                        items[i] = new StsColorListItem(StsColor.GREEN, vspName);
                    }
                }
                int[] selected = panel.getSelectedLeftIndices();

                // check to see if any of the currently selected have been removed
                if( selected != null )
                {
                	int nSelected = selected.length;
                	for( int i=0; i<selected.length; i++ )
                    	if( selected[i] >= nVsps ) nSelected--;
                    if( nSelected == 0 ) selected = null;
                    else
                    {
                    	int[] newSelected = new int[nSelected];
                        int iSelected = 0;
                        for( int i=0; i<nSelected; i++ )
	                    	if( selected[i] < nVsps )
                            {
                            	newSelected[iSelected] = selected[i];
                                iSelected++;
                            }
                        selected = new int[nSelected];
                        System.arraycopy(newSelected, 0, selected, 0, nSelected);
                    }

                }
                panel.setLeftItems(items);
                if( selected != null ) panel.setSelectedLeftIndices(selected);
            }
        }
    }

    private StsWell findWell(StsVsp vsp, int maxDistance)
    {
        StsClass wells = model.getCreateStsClass(StsWell.class);
        if(wells == null)
            return null;
        for(int i=0; i<wells.getSize(); i++)
        {
           StsWell well = (StsWell)wells.getElement(i);
           StsPoint pt1 = new StsPoint(well.getXOrigin(), well.getYOrigin(), 0.0f);
           StsPoint pt2 = new StsPoint(vsp.getXOrigin(), vsp.getYOrigin(), 0.0f);
           if(StsMath.distance(pt1, pt2) < maxDistance)
               return well;
        }
        return null;
    }

    private void initWellList()
    {
        int nWells = 0;
		StsClass wells = model.getCreateStsClass(StsWell.class);
        if(wells != null)
            nWells = wells.getSize();

        StsColorListItem[] items = new StsColorListItem[nWells+1];
        items[0] = new StsColorListItem(StsColor.GREEN, "New Well...");

        for( int i=1; i<nWells+1; i++ )
        {
            StsWell well = (StsWell)wells.getElement(i-1);
            items[i] = new StsColorListItem(well.getStsColor(), well.getName());
        }
        panel.setRightItems(items);
    }

    public boolean start()
    {
        enableFinish();
		refreshLists();

        // select all of the vsps
        if(vsps != null)
        {
        	int[] selected = new int[vsps.length];
            for( int i=0; i<vsps.length; i++ )
              selected[i] = i;
            panel.setSelectedLeftIndices(selected);
        }
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public void refreshLists()
    {
    	initWellList();
        initVspList();
        refreshButtons();
    }

    public String trimName(String fullString)
    {
    	int index = fullString.indexOf("(");
        if( index > 0 )
			return fullString.substring(0, index-1);
        return null;
    }

    public StsWell[] getSelectedWells()
    {
        StsWell[] wells = null;

        try
        {
            Object[] items = panel.getSelectedRightItems();
            int nItems = items == null ? 0 : items.length;
            if( nItems == 0) return null;

            wells = new StsWell[nItems];
            for( int i=0; i<nItems; i++ )
            {
                StsColorListItem item = (StsColorListItem)items[i];
                wells[i] = (StsWell)model.getObjectWithName(StsWell.class, item.getName());
            }
    	    return wells;
        }
        catch(Exception e)
        {
            StsException.outputException("StsVspWellAssign.getSelectedWells() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public void valueChanged(ListSelectionEvent e)
    {
    	refreshButtons();
    }

    public void refreshButtons()
    {
		Object[] vspItems = panel.getSelectedLeftItems();
		Object[] wellItems = panel.getSelectedRightItems();
    	int nVspItems = vspItems == null ? 0 : vspItems.length;
    	int nWellItems = wellItems == null ? 0 : wellItems.length;
        if((nVspItems > 0) && ((nVspItems == nWellItems) || (nWellItems == 1)))
           panel.enableLeftButton(true);

		if(nVspItems > 0)
        {
			for(int i=0; i<nVspItems; i++)
            {
            	StsColorListItem item = (StsColorListItem) vspItems[i];
            	StsVsp vsp = null;
                try { vsp = (StsVsp)model.getObjectWithName(StsVsp.class, trimName(item.getName())); }
                catch(Exception ex) { }
                if( vsp != null && vsp.getWell() != null )
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
			Object[] vspItems = panel.getSelectedLeftItems();
			Object[] wellItems = panel.getSelectedRightItems();
	    	int nVspItems = vspItems == null ? 0 : vspItems.length;
    		int nWellItems = wellItems == null ? 0 : wellItems.length;

            // Dont allow New Well and other wells to be selected
            if(nWellItems > 1)
            {
                for (int i = 0; i < nWellItems; i++)
                {
                    if (((StsColorListItem) wellItems[i]).getName().equals("New Well..."))
                    {
                        new StsMessage(wizard.frame, StsMessage.WARNING,
                            "Cannot select New Well with any other exisitng wells");
                        panel.clearRightSelections();
                        return;
                    }
                }
            }

            // If the New Well item was selected
            if((nWellItems == 1) && (((StsColorListItem)wellItems[0]).getName().equals("New Well...")))
            {
                System.out.println("Define new well...");
                // Show dialog to add a well
                WellDefineDialog wDialog = new WellDefineDialog(wizard.frame);
                wDialog.setVisible(true);
                StsWell well = wDialog.getWell();
				if (well == null) return ; // jbw
                // Add new well to right list
                refreshLists();

                Object[] wells = panel.getRightItems();
                for(int i=0; i<wells.length; i++)
                {
                    if(((StsColorListItem)wells[i]).getName().equals(well.getName()))
                      panel.setSelectedRightIndices(new int[] {i});
                }
                // Set new well as only selected well
                wellItems =  panel.getSelectedRightItems();
                if(wellItems.length == 0)
                    return;
            }

            // If multiple select on both list, must be same number
            if(nVspItems == nWellItems )
            {
                for( int i=0; i<nVspItems; i++ )
                {
                    StsColorListItem item = (StsColorListItem)vspItems[i];
                    StsVsp vsp = null;
                    try { vsp = (StsVsp)model.getObjectWithName(StsVsp.class, trimName(item.getName())); }
                    catch(Exception ex)
                    {
                        new StsException(StsException.WARNING, "Unable to find well" + ex);
                    }

                    item = (StsColorListItem) wellItems[i];
                    StsWell well = (StsWell)model.getObjectWithName(StsWell.class, item.getName());

                    if( vsp != null && well != null )
                    {
                        setVspDisplayFlags(well);
                        vsp.setWell(well);
                    }
                }
                refreshLists();
            }

            // If mulitple VSPs and one well
            else if(nWellItems == 1)
            {
                StsColorListItem wellItem = (StsColorListItem) wellItems[0];
                StsWell well = (StsWell)model.getObjectWithName(StsWell.class, wellItem.getName());

                for( int i=0; i<nVspItems; i++ )
                {
                    StsColorListItem item = (StsColorListItem) vspItems[i];
                    StsVsp vsp = null;
                    try { vsp = (StsVsp)model.getObjectWithName(StsVsp.class, trimName(item.getName())); }
                    catch(Exception ex) { }

                    if( vsp != null && well != null )
                    {
                        setVspDisplayFlags(well);
                        vsp.setWell(well);
                    }
                }
                refreshLists();
            }
        }

        else if( e.getActionCommand().equals(">") )
        {
			Object[] vspItems = panel.getSelectedLeftItems();
	    	int nVspItems = vspItems == null ? 0 : vspItems.length;
            for( int i=0; i<nVspItems; i++ )
            {
                StsColorListItem item = (StsColorListItem) vspItems[i];
                StsVsp vsp = null;
                try { vsp = (StsVsp)model.getObjectWithName(StsVsp.class, trimName(item.getName())); }
                catch(Exception ex) { }
                if( vsp != null )
                    vsp.setWell(null);
            }
            refreshLists();
        }
    }

    public void setVspDisplayFlags(StsWell well)
    {
        // Determine if any other vsps are assigned to this well?
        StsObject[] vsps = (StsObject[]) model.getObjectList(StsVsp.class);
        for(int i=0; i<vsps.length; i++)
        {
            if (((StsVsp)vsps[i]).getWell() == well)
                ((StsVsp)vsps[i]).setDisplayVSPs(false);
        }
    }

    class WellDefineDialog extends JDialog
    {
        StsWell well = null;

        StsJPanel panel = StsJPanel.addInsets();
        GridBagLayout gridBagLayout1 = new GridBagLayout();

        Frame frame = null;
        StsStringFieldBean nameBean = null;
        String wellName = "wellName";

        StsGroupBox defineTopBox = new StsGroupBox("Define Top Location");
        StsFloatFieldBean topXBean = null;
        StsFloatFieldBean topYBean = null;
        StsFloatFieldBean topZBean = null;
        StsFloatFieldBean topTBean = null;
        float topX = 0.0f;
        float topY = 0.0f;
        float topZ = 0.0f;
        float topT = 0.0f;

        StsGroupBox defineBtmBox = new StsGroupBox("Define Bottom Location");
        StsFloatFieldBean btmXBean = null;
        StsFloatFieldBean btmYBean = null;
        StsFloatFieldBean btmZBean = null;
        StsFloatFieldBean btmTBean = null;
        float btmX = 0.0f;
        float btmY = 0.0f;
        float btmZ = 0.0f;
        float btmT = 0.0f;

        StsButton okBtn = new StsButton();
        StsButton cancelBtn = new StsButton();

        public WellDefineDialog(Frame frame)
        {
            super(frame, "Well Definition", true);
            this.setLocation(frame.getLocation());
            this.frame = frame;
            try
            {
                constructBeans();
                jbInit();
                pack();
            }
            catch(Exception ex)
            {
                ex.printStackTrace();
            }
        }

        private void constructBeans()
        {
            nameBean = new StsStringFieldBean(this, "wellName", true, "Well Name:");
            topXBean = new StsFloatFieldBean(this, "topX", true, "X:", false);
            topYBean = new StsFloatFieldBean(this, "topY", true, "Y:", false);
            topZBean = new StsFloatFieldBean(this, "topZ", true, "Z:", false);
            topTBean = new StsFloatFieldBean(this, "topT", true, "Time:", false);

            btmXBean = new StsFloatFieldBean(this, "btmX", true, "X:", false);
            btmYBean = new StsFloatFieldBean(this, "btmY", true, "Y:", false);
            btmZBean = new StsFloatFieldBean(this, "btmZ", true, "Z:", false);
            btmTBean = new StsFloatFieldBean(this, "btmT", true, "Time:", false);

            okBtn = new StsButton("Ok", "Create well using input values", this, "okProcess", null);
            cancelBtn = new StsButton("Cancel", "Cancel well creation", this, "cancelProcess", null);

        }

        private void jbInit() throws Exception
        {
            panel.gbc.fill = panel.gbc.HORIZONTAL;
            panel.setBorder(BorderFactory.createEtchedBorder());
            panel.setLayout(gridBagLayout1);

            panel.add(nameBean);

            panel.gbc.gridwidth = 2;
            defineTopBox.addEndRow(topXBean);
            defineTopBox.addEndRow(topYBean);
            defineTopBox.addEndRow(topZBean);
            defineTopBox.addEndRow(topTBean);
            panel.add(defineTopBox);

            defineBtmBox.addEndRow(btmXBean);
            defineBtmBox.addEndRow(btmYBean);
            defineBtmBox.addEndRow(btmZBean);
            defineBtmBox.addEndRow(btmTBean);
            panel.add(defineBtmBox);

            panel.gbc.gridwidth = 1;
            panel.addToRow(okBtn);
            panel.addEndRow(cancelBtn);

            this.getContentPane().add(panel);
        }

        public String getWellName() { return wellName; }
        public float getTopX() { return topX; }
        public float getTopY() { return topY; }
        public float getTopZ() { return topZ; }
        public float getTopT() { return topT; }

        public float getBtmX() { return btmX; }
        public float getBtmY() { return btmY; }
        public float getBtmZ() { return btmZ; }
        public float getBtmT() { return btmT; }

        public void setWellName(String value) { wellName = value; }
        public void setTopX(float value) { topX = value; }
        public void setTopY(float value) { topY = value; }
        public void setTopZ(float value) { topZ = value; }
        public void setTopT(float value) { topT = value; }

        public void setBtmX(float value) { btmX = value; }
        public void setBtmY(float value) { btmY = value; }
        public void setBtmZ(float value) { btmZ = value; }
        public void setBtmT(float value) { btmT = value; }

        public void okProcess()
        {
            well = createWell();
            this.setVisible(false);
        }

        public void cancelProcess()
        {
            well = null;
            this.hide();
        }

        public StsWell getWell() { return well; }

        public StsWell createWell()
        {
            StsLogVector[] vectors = new StsLogVector[5];
            StsLogCurve tdCurve = null;

            model.disableDisplay();
            StsWell well = new StsWell(wellName, false);

            vectors[0] = new StsLogVector(StsLogVector.X, new float[] {getTopX(), getBtmX()} );
            vectors[0].setUnits(StsParameters.DIST_FEET);
            vectors[1] = new StsLogVector(StsLogVector.Y, new float[] {getTopY(), getBtmY()} );
            vectors[1].setUnits(StsParameters.DIST_FEET);
            vectors[2] = new StsLogVector(StsLogVector.DEPTH, new float[] {getTopZ(), getBtmZ()} );
            vectors[2].setUnits(StsParameters.DIST_FEET);
            vectors[3] = new StsLogVector(StsLogVector.MDEPTH, new float[] {getTopZ(), getBtmZ()} );
            vectors[3].setUnits(StsParameters.DIST_FEET);
            vectors[4] = new StsLogVector(StsLogVector.TIME, new float[] {getTopT(), getBtmT()} );
            vectors[4].setUnits(StsParameters.TIME_MSECOND);
            tdCurve = new StsLogCurve(vectors[3], vectors[2], vectors[4], 0);

            well.constructWellDevCurves(vectors, StsParameters.nullValue, tdCurve);
            well.addLogCurve(tdCurve);
            well.computePoints();
            well.setXOrigin(getTopX());
            well.setYOrigin(getTopY());
            well.setZDomainSupported(well.TD_TIME_DEPTH);
            well.addToModel();
            well.addToProject();

            StsProject project = wizard.getModel().getProject();
            project.adjustBoundingBoxes(true, true);
            project.checkAddUnrotatedClass(StsWell.class);
            project.rangeChanged();
            StsView[] views = model.win3d.getDisplayedViews();
            for(StsView view : views)
                view.glPanel3d.setDefaultView();
            model.enableDisplay();
            model.win3dDisplay();

            return well;
        }
    }
}
