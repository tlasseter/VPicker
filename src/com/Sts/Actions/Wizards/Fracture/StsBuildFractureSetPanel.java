package com.Sts.Actions.Wizards.Fracture;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsBuildFractureSetPanel extends StsJPanel
{
    private StsWizard wizard;
    private StsBuildFractureSet wizardStep;
    public StsGroupBox addBox = new StsGroupBox("Add Edges and Sections");
    public StsButton endEdgeButton;
    public StsButton endSectionButton;
    public StsGroupBox deleteBox  = new StsGroupBox("Delete Edges");
    public StsButton deleteCurrentEdgeButton;
    public StsButton deleteSectionButton;

    public StsBuildFractureSetPanel(StsWizard wizard, StsBuildFractureSet wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
        	buildPanel();
        }
        catch(Exception e) {}
    }

    public void initialize()
    {
    }

    void buildPanel() throws Exception
    {
        endEdgeButton = new StsButton("End Edge", "End picks on current edge.", wizardStep, "endEdge");
        endSectionButton = new StsButton("End Section", "End current edge and complete section.", wizardStep, "endSection");
        deleteCurrentEdgeButton = new StsButton("Delete Edge", "Delete current edge.", wizardStep, "deleteCurrentEdge");
        //deleteSectionButton = new StsButton("Delete Section", "Delete current section.", wizardStep, "deleteSection");
        
        addBox.addToRow(endEdgeButton);
        addBox.addEndRow(endSectionButton);
        gbc.fill = gbc.HORIZONTAL;
        add(addBox);
        deleteBox.addToRow(deleteCurrentEdgeButton);
        //deleteBox.addToRow(deleteSectionButton);        
        add(deleteBox);
    }
}
