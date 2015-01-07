
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.AncillaryData;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

public class StsAncillaryAssign extends StsWizardStep implements ActionListener, ListSelectionListener
{
    StsHeaderPanel header = new StsHeaderPanel();
	JPanel mainPanel = new JPanel(new BorderLayout());
    StsAncillaryData[] docs = null;
	StsDoubleListPanel panel = new StsDoubleListPanel();
    StsAncillaryDataWizard wizard = null;

    public StsAncillaryAssign(StsAncillaryDataWizard wizard)
    {
        super(wizard);
        this.wizard = wizard;
        setPanels(mainPanel, header);
        mainPanel.add(panel);
        panel.setLeftListRenderer(new StsColorListRenderer());
        panel.setRightListRenderer(new StsColorListRenderer());
        panel.setLeftTitle("Selected Files:");
        panel.setRightTitle("Available Wells:");
        panel.addActionListener(this);
        panel.addListSelectionListener(this);
        panel.getLeftButton().setToolTipText("Assign well to selected File");
        panel.getRightButton().setToolTipText("Remove well assignment from selected File");
        panel.setPreferredSize(new Dimension(400,300));

        header = (StsHeaderPanel) getHdrContainer();

        header.setTitle("Ancillary Data Load");
        header.setSubtitle("Assign Wells to Ancillary Data");
        header.setLink("http://www.s2ssystems.com/Marketing/AppLinks.html#AncillaryData");                
        header.setInfoText(wizardDialog,"(1) Select a file from left list.\n" +
                                  "(2) Select the related well from the right list.\n" +
                                  "   ***** Well list is generated from all wells loaded in project *****\n" +
                                  "   ***** Psuedo wells can be created prior to assignment *****\n" +
                                  "(3) Press New Well.. Button to create vertical psuedo well\n" +
                                  "(4) Press the < button to assign the well to the selected file or\n" +
                                  "    Press the > button to unassign the well from the selected file.\n" +
                                  "(5) Press the Next>> Button to complete the ancillary data load.");
    }

    public void initDocumentList()
    {
        docs = wizard.getAncillaryDataObjects();

        if(docs != null )
        {
            int nDocs = docs == null ? 0 : docs.length;
            if( nDocs > 0 )
            {
                StsColorListItem[] items = new StsColorListItem[nDocs];
                for( int i=0; i<nDocs; i++ )
                {
                    StsAncillaryData doc = (StsAncillaryData)docs[i];
                    String docName = doc.getName();
                    StsWell well = doc.getWell();
                    if(well == null)
                        doc.setWell(well);
                    if (well == null)
                    {
                        docName += " ( )";
                        items[i] = new StsColorListItem(doc);
                    }
                    else
                    {
                        docName += " (" + well.getName() + ")";
                        items[i] = new StsColorListItem(doc);
                    }
                    items[i].setName(docName);
                }
                int[] selected = panel.getSelectedLeftIndices();

                // check to see if any of the currently selected have been removed
                if( selected != null )
                {
                	int nSelected = selected.length;
                	for( int i=0; i<selected.length; i++ )
                    	if( selected[i] >= nDocs ) nSelected--;
                    if( nSelected == 0 ) selected = null;
                    else
                    {
                    	int[] newSelected = new int[nSelected];
                        int iSelected = 0;
                        for( int i=0; i<nSelected; i++ )
	                    	if( selected[i] < nDocs )
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

        // select all of the docs
        if(docs != null)
        {
        	int[] selected = new int[docs.length];
            for( int i=0; i<docs.length; i++ )
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
        initDocumentList();
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
            StsException.outputException("StsAncillaryAssign.getSelectedWells() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public void valueChanged(ListSelectionEvent e)
    {
    	refreshButtons();
    }

    public void refreshButtons()
    {
		Object[] docItems = panel.getSelectedLeftItems();
		Object[] wellItems = panel.getSelectedRightItems();
    	int nDocItems = docItems == null ? 0 : docItems.length;
    	int nWellItems = wellItems == null ? 0 : wellItems.length;
        if((nDocItems > 0) && ((nDocItems == nWellItems) || (nWellItems == 1)))
           panel.enableLeftButton(true);

		if(nDocItems > 0)
        {
			for(int i=0; i<nDocItems; i++)
            {
            	StsColorListItem item = (StsColorListItem) docItems[i];
                StsAncillaryData doc = (StsAncillaryData)item.getObject();
                if( doc != null && doc.getWell() != null )
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
			Object[] docItems = panel.getSelectedLeftItems();
			Object[] wellItems = panel.getSelectedRightItems();
	    	int nDocItems = docItems == null ? 0 : docItems.length;
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
                wDialog.show();
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
            if(nDocItems == nWellItems )
            {
                for( int i=0; i<nDocItems; i++ )
                {
                    StsColorListItem item = (StsColorListItem)docItems[i];
                    StsAncillaryData doc = (StsAncillaryData)item.getObject();
                    item = (StsColorListItem) wellItems[i];
                    StsWell well = (StsWell)model.getObjectWithName(StsWell.class, item.getName());

                    if( doc != null && well != null )
                        doc.setWell(well);
                }
                refreshLists();
            }

            // If mulitple Files and one well
            else if(nWellItems == 1)
            {
                StsColorListItem wellItem = (StsColorListItem) wellItems[0];
                StsWell well = (StsWell)model.getObjectWithName(StsWell.class, wellItem.getName());

                for( int i=0; i<nDocItems; i++ )
                {
                    StsColorListItem item = (StsColorListItem) docItems[i];
                    StsAncillaryData doc = (StsAncillaryData)item.getObject();
                    if( doc != null && well != null )
                        doc.setWell(well);
                }
                refreshLists();
            }
        }

        else if( e.getActionCommand().equals(">") )
        {
			Object[] docItems = panel.getSelectedLeftItems();
	    	int nDocItems = docItems == null ? 0 : docItems.length;
            for( int i=0; i<nDocItems; i++ )
            {
                StsColorListItem item = (StsColorListItem) docItems[i];
                StsAncillaryData doc = (StsAncillaryData)item.getObject();
                if( doc != null )
                    doc.setWell(null);
            }
            refreshLists();
        }
    }

    class WellDefineDialog extends JDialog
    {
        StsWell well = null;

        StsJPanel panel = new StsJPanel();
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
            this.hide();
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
//            StsGLPanel3d glPanel3d = (StsGLPanel3d)wizard.getGlPanel();
//            glPanel3d.cursor3d.initialize();
//            model.win3d.cursor3dPanel.setSliderValues();
//            model.glPanel3d.setDefaultView();
            model.enableDisplay();
            model.win3dDisplay();

            return well;
        }
    }
}
