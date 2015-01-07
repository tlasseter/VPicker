
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001i
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Types.*;

import java.awt.*;

public class StsBuildSectionsFromEdges extends StsAction implements Runnable
{
    boolean success = false;
    boolean debug;

    static final int MINUS = StsParameters.MINUS;
    static final int PLUS = StsParameters.PLUS;

    static final int FAULT = StsParameters.FAULT;
    static final int AUXILIARY = StsParameters.AUXILIARY;

 	public StsBuildSectionsFromEdges(StsActionManager actionManager)
    {
        super(actionManager);
        statusArea.setTitle("Build:");
        logMessage("Making fault sections from surface cuts...");
        debug = model.getBooleanProperty("debugFrame");
    }

    public boolean start()
    {
        return true;
    }

	public void run()
    {
        try
        {
			buildFrame();
			actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            StsException.outputException("StsBuildSectionsFromEdges.run() failed.", e, StsException.WARNING);
        }
    }

    public boolean end()
    {
    	statusArea.textOnly();
    	if( success )
           	logMessage("Fault section(s) built successfully.");
        else
        	logMessage("Error building fault section(s).");
        return success;
    }

	private void buildFrame()
    {
        StsCursor cursor = new StsCursor(model.win3d, Cursor.WAIT_CURSOR);

        if(!constructSectionsFromEdges())
        {
            success = false;
            cursor.restoreCursor();
            return;
        }

        intersectSections();
        success = true;
        cursor.restoreCursor();
        model.win3d.win3dDisplay();
    }

    private boolean constructSectionsFromEdges()
    {
        StsSectionEdge edge;
        StsSurfaceVertex vertex;
        StsLine firstLine, lastLine;
        StsSection section = null;

        try
        {
            StsClass edges = model.getCreateStsClass(StsSurfaceEdge.class);
            int nEdges = edges.getSize();

            int nSectionsBuilt = 0;
		    boolean constructionOK = true;
            for(int n = 0; n < nEdges; n++)
            {
                edge = (StsSectionEdge)edges.getElement(n);
                if( edge.getType() == StsSectionEdge.REFERENCE ) continue;
                if(edge.getSection() != null) continue;
                float dZMax = 2.0f * model.getProject().getZorTInc();

                boolean isCurved = edge.isCurved();

                vertex = edge.getPrevVertex();
                firstLine = vertex.getSectionLine();
                if(firstLine == null)
                {
                    firstLine = StsLine.buildVertical(vertex, StsParameters.FAULT);
                    firstLine.insertLineVertex(vertex, false, dZMax, true);  // share w/ edge
  //                  firstLine.insertSurfaceVertex(vertex);
                }
                vertex = edge.getNextVertex();
                lastLine = vertex.getSectionLine();
                if(lastLine == null)
                {
                    lastLine = StsLine.buildVertical(vertex, StsParameters.FAULT);
                    lastLine.insertLineVertex(vertex, false, dZMax, true); // share w/ edge
  //                  lastLine.insertSurfaceVertex(vertex);
                }
                section = StsSection.constructor(edge, edge.getType(), firstLine, lastLine, StsSectionGeometry.CURVED);
				if(section == null)
				{
					logMessage("Build Sections failed for section: " + section.getLabel() + " See error log for details.");
					constructionOK = false;
				}
                nSectionsBuilt++;
            }
            return constructionOK;
        }
        catch(Exception e)
        {
            StsException.outputException("StsBuildSectionsFromEdges.constructSectionsFromEdges() failed.",
				e, StsException.WARNING);
            return false;
        }
    }

    private void intersectSections()
    {
        StsSectionEdge edge;

        StsClass edges = model.getCreateStsClass(StsSurfaceEdge.class);
        int nEdges = edges.getSize();

        try
        {
            for(int n = 0; n < nEdges; n++)
            {
                edge = (StsSectionEdge)edges.getElement(n);
                byte type = edge.getType();
                if(type != FAULT && type != AUXILIARY) continue;
                intersectSections(edge, MINUS);
                intersectSections(edge, PLUS);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsBuildSectionsFromEdges.intersectSections() failed.",
                e, StsException.WARNING);
            return;
        }
    }

    private boolean intersectSections(StsSectionEdge edge, int end)
    {
        StsSection section;
        StsLine line;
        StsSurfaceVertex vertex;

        try
        {
            section = edge.getSection();
            if(section == null) return false;

            if(end == MINUS)
            {
                vertex = edge.getPrevVertex();
                line = section.getFirstLine();
            }
            else
            {
                vertex = edge.getNextVertex();
                line = section.getLastLine();
            }

            StsSectionEdge connectEdge = vertex.getSectionEdge();
            if(connectEdge == null) return false;

            StsSection connectSection = connectEdge.getSection();

            line.addOnSection(connectSection);

            section.addLineToSectionSide(end);

            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsBuildSectionsFromEdges.intersectSections() failed.",
                                         e, StsException.WARNING);
            return false;
        }
    }
}
