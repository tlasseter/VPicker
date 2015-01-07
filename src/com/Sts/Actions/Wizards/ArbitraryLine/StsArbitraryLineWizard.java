
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.ArbitraryLine;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsArbitraryLineWizard extends StsWizard
{
    private ArrayList lines = new ArrayList();

    public StsArbitraryLineSelect lineSelect = null;
    public StsArbitraryLineProperty lineProperty = null;
    public StsArbitraryLineIntersect lineIntersect = null;
    public StsArbitraryLineDefine lineDefine = null;

    private StsArbitraryLine selectedLine = null;

    private StsWizardStep[] mySteps =
    {
        lineSelect = new StsArbitraryLineSelect(this),
        lineProperty = new StsArbitraryLineProperty(this),
        lineDefine = new StsArbitraryLineDefine(this),
        lineIntersect = new StsArbitraryLineIntersect(this),
    };


    private byte pickMode = NO_PICK_MODE;

    static final byte NO_PICK_MODE = 0;
    static final byte LINE_CREATE = 1;
    static final byte LINE_MOVE = 2;
    static final byte LINE_DELETE = 3;
    static final byte LINE_ROTATE = 4;
    static final byte INTERSECT_ADD = 5;
    static final byte INTERSECT_DELETE = 6;

    public StsArbitraryLineWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Arbitrary Line Setup");
        dialog.getContentPane().setSize(300, 400);

        if(!super.start()) return false;

        StsArbitraryLineClass lineClass = StsArbitraryLine.getArbitraryLineClass();
        int nLines = lineClass.getSize();

        StsArbitraryLine currentLine = lineClass.getCurrentArbitraryLine();
        if(currentLine == null && nLines > 0)
        {
            currentLine = (StsArbitraryLine)lineClass.getLast();
            //model.addSelectedObject(currentLine);
            lineClass.setCurrentObject(currentLine);
        }
        return true;
    }

    public void createNewLine()
    {
        gotoStep(lineProperty);
    }

    public void setSelectHinge(boolean val)
    {
        setPickMode(LINE_CREATE);
    }

    public StsObjectRefList getHinges()
    {
        return null;
    }

    public StsObjectRefList getIntersections()
    {
        return null;
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
       gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == lineSelect)
        {
            selectedLine = lineSelect.getSelectedLine();
            if(selectedLine == null)
            {
                disableFinish();
                new StsMessage(frame, StsMessage.INFO, "No line selected.");
                return;
            }
            else
            {
                enableFinish();
                model.setCurrentObject(selectedLine);
                super.finish();
            }
        }
        else
        {
            gotoNextStep();
        }

    }

    public void finish()
    {
        super.finish();
    }

    public void setPickMode(byte pickMode) { this.pickMode = pickMode; }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        switch(pickMode)
        {
            case NO_PICK_MODE:
                break;
            case LINE_CREATE:
                return createLine(mouse, glPanel);
            case LINE_MOVE:
                break;
            case LINE_DELETE:
                break;
            case LINE_ROTATE:
                break;
            case INTERSECT_ADD:
                break;
             case INTERSECT_DELETE:
        }
        model.win3dDisplay();
        return true;
    }

    private boolean createLine(StsMouse mouse, StsGLPanel glPanel)
    {
        StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getStsClass(StsSurface.class);
        StsSurface surface = surfaceClass.getTopVisibleSurface();
        if(surface != null)
        {
            int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

            if (leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
            {
                surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
            }
            else // mouse.RELEASED: add a vertical line at this point
            {
                StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
                if(gridPoint == null) return true;
                StsLine line = StsLine.buildVertical(gridPoint, StsParameters.FAULT);
                if(line != null)
                {
                    lines.add(line);
                    int nLines = lines.size();
                    if(nLines > 1)
                    {
                        StsLine prevLine = (StsLine)lines.get(nLines-2);
                        try
                        {
                            new StsSection(null, StsSection.AUXILIARY, prevLine, line, StsSectionGeometry.RIBBON);
                        }
                        catch(Exception e)
                        {
                            new StsMessage(model.win3d, StsMessage.WARNING, "Failed to construct section between " +
                                           prevLine.getLabel() + " and " + line.getLabel());
                        }
                    }
                    model.win3d.win3dDisplay();
                }
            }
            return true;
        }
        return false;
    }
}
