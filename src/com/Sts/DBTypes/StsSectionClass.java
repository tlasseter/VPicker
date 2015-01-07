package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

public class StsSectionClass extends StsClass implements StsSerializable, StsClassDisplayable, StsClassTimeDepthDisplayable
{
    protected boolean displaySections = true;

    protected boolean displayFaultSections = true;
    protected boolean displayAuxiliarySections = true;
    protected boolean displayBoundarySections = true;
    protected boolean displayFractureSections = true;
    
    transient protected boolean displayFaultSectionEdges = true;
    transient protected boolean displayAuxiliarySectionEdges = true;
    transient protected boolean displayBoundarySectionEdges = true;
    transient protected boolean displayFractureSectionEdges = true;
    
    /** zDomain(s) supported; set to two domains (time & depth or approxDepth & depth) if a velocity model is available and has been applied. */
     protected byte zDomainSupported = StsParameters.TD_NONE;

	static public final StsColor defaultFaultColor = StsColor.PURPLE;

    public StsSectionClass()
    {
    }

    public void initializeDisplayFields()
    {
//		initColors(StsSection.displayFields);
//		initColors(StsSection.faultDisplayFields);

        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(this, "displaySections", "Enable")
        };
    }

    public boolean getDisplaySectionEdges(int type)
    {
        switch (type)
        {
            case StsParameters.FAULT:
                return displayFaultSectionEdges;
            case StsParameters.AUXILIARY:
                return displayAuxiliarySectionEdges;
            case StsParameters.BOUNDARY:
                return displayBoundarySectionEdges;
            case StsParameters.FRACTURE:
                return displayFractureSectionEdges;                
            default:
                return false;
        }
    }

    public void setDisplaySections(boolean b)
    {
        displaySections = b;
        displayFaultSections = b;
        displayAuxiliarySections = b;
        displayBoundarySections = b;
    }

    public void setDisplaySectionEdges(boolean b)
    {
        displayFaultSectionEdges = b;
        displayAuxiliarySectionEdges = b;
        displayBoundarySectionEdges = b;
        setDrawEdges(b);
    }

    public void setDrawEdges(boolean b)
    {
        for(int n = 0; n < getSize(); n++)
        {
            StsSection section = (StsSection)getElement(n);
            section.drawEdges = b;
        }
    }

    public boolean getDisplaySections() { return displaySections; }

    public boolean getDisplayFaultSections() { return displayFaultSections; }
    public boolean getDisplayBoundarySections() { return displayBoundarySections; }
    public boolean getDisplayAuxiliarySections() { return displayAuxiliarySections; }
    public boolean getDisplayFractureSections() { return displayFractureSections; }
    
    public void setDisplayFaultSections(boolean b) { displayFaultSections = b; }
    public void setDisplayBoundarySections(boolean b) { displayBoundarySections = b; }
    public void setDisplayAuxiliarySections(boolean b) { displayAuxiliarySections = b; }
    public void setDisplayFractureSections(boolean b) { displayFractureSections = b; }
    
    public boolean getDisplayFaultSectionEdges() { return displayFaultSectionEdges; }
    public boolean getDisplayBoundarySectionEdges() { return displayBoundarySectionEdges; }
    public boolean getDisplayAuxiliarySectionEdges() { return displayAuxiliarySectionEdges; }
    public boolean getDisplayFractureSectionEdges() { return displayFractureSectionEdges; }
    
    public void setDisplayFaultSectionEdges(boolean b) { displayFaultSectionEdges = b; }
    public void setDisplayBoundarySectionEdges(boolean b) { displayBoundarySectionEdges = b; }
    public void setDisplayAuxiliarySectionEdges(boolean b) { displayAuxiliarySectionEdges = b; }
    public void setDisplayFractureSectionEdges(boolean b) { displayFractureSectionEdges = b; }
    
    public void displayClass(StsGLPanel3d glPanel3d)
    {
        if(!displayBoundarySections && !displayFaultSections && !displayAuxiliarySections && !displayFractureSections) return;


        int nSections = getSize();
        for(int n = 0; n < nSections; n++)
        {
            StsSection section = (StsSection)getElement(n);
            switch(section.getType())
            {
                case StsSection.BOUNDARY:
                    if(displayBoundarySections) section.display(glPanel3d, displayBoundarySectionEdges);
                    break;
                case StsSection.FAULT:
                    if(displayFaultSections) section.display(glPanel3d, displayFaultSectionEdges);
                    break;
                case StsSection.FRACTURE:
                    if(displayFractureSections) section.display(glPanel3d, displayFractureSectionEdges);
                    break;                    
                case StsSection.AUXILIARY:
                    if(displayAuxiliarySections) section.display(glPanel3d, displayAuxiliarySectionEdges);
            }
        }
    }
    public StsObject[] getVisibleSectionEdges()
    {
        StsSectionEdge[] edges = new StsSectionEdge[0];

        int nSections = getSize();
        for (int n = 0; n < nSections; n++)
        {
            StsSection section = (StsSection)getElement(n);

            if (!section.getIsVisible()) continue;
            if (section.getType() == StsSection.BOUNDARY) continue;

            StsObjectRefList sectionEdges = section.getSectionEdges();
            int nEdges = sectionEdges.getSize();
            for (int e = 0; e < nEdges; e++)
            {
                StsSectionEdge edge = (StsSectionEdge)sectionEdges.getElement(e);
                StsXYSurfaceGridable surface = edge.getSurface();
                if (surface != null && !surface.getIsVisible()) continue;
                edges = (StsSectionEdge[])StsMath.arrayAddElement(edges, edge);
            }
        }
        return edges;
    }

    public void checkSections()
    {
        int nSections = getSize();

        for (int n = nSections - 1; n >= 0; n--)
        {
            StsSection section = (StsSection) getElement(n);
            section.checkDelete();
        }
    }

    /** Initialize any sections for which both defining lines are initialized */
    public void initSectionsNotOnSections()
    {
        int nSections = getSize();

        for (int n = 0; n < nSections; n++)
        {
            StsSection section = (StsSection) getElement(n);
            if (section != null && !section.initialized)
            {
                StsLine[] lines = section.getLineList();
                if (lines[0].initialized && lines[1].initialized)
                {
                    section.initialize();
                }
            }
        }
    }

    public boolean initSections()
    {
        boolean allInitialized = true;

        int nSections = getSize();

        for (int n = 0; n < nSections; n++)
        {
            StsSection section = (StsSection) getElement(n);
            if (section != null && !section.initialized)
            {
                if (!section.initialize())
                {
                    allInitialized = false;
                }
            }
        }

        return allInitialized;
    }

    /** if there are no fault or auxiliary sections, set the project.modelZDomain to NONE; otherwise leave as is.
     *  If none, subsequently when a fault or aux section is built, the modelZDomain will be initialized
     *  to the then current zDomain.
     */
    public void checkSetModelZDomain()
    {
        for(int n = 0; n < getSize(); n++)
        {
            StsSection section = (StsSection)getElement(n);
            if(section.getType() != StsSection.BOUNDARY) return;
        }
        currentModel.getProject().setModelZDomain(StsProject.TD_NONE);
    }
}
/*
    static public final StsFieldBean[] boundaryFieldBeans = {
        new StsBooleanFieldBean(StsModel.class, DISPLAY_BOUNDARY, "Enable"),
        new StsBooleanFieldBean(StsModel.class,
                                   DISPLAY_BOUNDARY_SECTION_EDGES, "Edges")
    };
    static public final StsFieldBean[] faultFieldBeans = {
        new StsBooleanFieldBean(StsModel.class, DISPLAY_FAULTS, "Enable"),
        new StsBooleanFieldBean(StsModel.class, DISPLAY_FAULT_SECTION_EDGES,
                                   "Edges")
    };
    static public final StsFieldBean[] auxiliaryFieldBeans = {
        new StsBooleanFieldBean(StsModel.class, DISPLAY_AUXILIARIES,
                                   "Enable"),
        new StsBooleanFieldBean(StsModel.class,
                                   DISPLAY_AUXILIARY_SECTION_EDGES, "Edges")
    };
*/
