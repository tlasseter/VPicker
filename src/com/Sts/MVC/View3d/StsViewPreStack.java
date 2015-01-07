package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Base class used to present two-dimensional data.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.Actions.Wizards.VelocityAnalysis.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Reflect.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;

abstract public class StsViewPreStack extends StsView2d implements StsSerializable
{
    /**
     * Current prestack volume being viewed
     */
    transient public StsPreStackLineSet lineSet = null;
    /**
     * Model prestack volume class which maintains a list of all crossplot objects
     */
    transient public StsPreStackLineSetClass lineSetClass = null;
    transient int pickedIndex = -1;
    transient StsPoint pickedPoint;
    /**
     * Has the prestack volume changed and therefore needs to be re-drawn
     */
    transient public boolean volumeChanged = true;

    //    public int currentRow = -1;
    //    public int currentCol = -1;
    /**
     * velocity profile currently being drawn
     */
    protected transient StsVelocityProfile currentVelocityProfile;
    transient public StsSuperGather superGather;
    static final byte nullByte = StsParameters.nullByte;
    static final boolean debug = false;

    public StsViewPreStack()
    {
    }

    public StsViewPreStack(StsGLPanel3d glPanel3d)
    {
        super(glPanel3d);
    }

    public void setDefaultView(byte axis)
    {
        //super.setDefaultView();
        axisRanges = StsMath.copyFloatArray(totalAxisRanges);
        float tpi = lineSet.getStsPreStackSeismicClass().defaultWiggleDisplayProperties.getTracesPerInch();
        float ips = lineSet.getStsPreStackSeismicClass().defaultWiggleDisplayProperties.getInchesPerSecond();
        lineSet.getWiggleDisplayProperties().setTracesPerInch(tpi);
        lineSet.getWiggleDisplayProperties().setInchesPerSecond(ips);
        //        recomputeZoomedScaling(axisRanges, axis);
        computePixelScaling();
        repaint();
    }

    public void setDefaultView()
    {
        setDefaultView(BOTH);
    }

    public void setLineSet(StsPreStackLineSet lineSet)
    {
        this.lineSet = lineSet;
        lineSetClass = (StsPreStackLineSetClass) model.getCreateStsClass(lineSet);
        superGather = lineSet.getSuperGather(glPanel3d.window);
    }

    public void pickOnSeismicView(StsGLPanel3d glPanel3d, StsVelocityProfile velocityProfile)
    {
        if (!velocityProfile.isPersistent()) return;
        GL gl = glPanel3d.getGL();
        gl.glLineWidth(4.0f);
        StsPoint[] displayedPoints = getDisplayedProfilePoints(velocityProfile);
        pickVertices(gl, displayedPoints);
        //        pickMutes(gl, velocityProfile);
        pickLine(gl, displayedPoints);
    }

    public StsPoint[] getDisplayedProfilePoints(StsVelocityProfile velocityProfile)
    {
        return velocityProfile.getProfilePoints();
    }

    public void pickVertices(GL gl, StsPoint[] displayedPoints)
    {
        gl.glInitNames();
        gl.glPushName(StsVelocityAnalysisEdit.TYPE_GATHER_VERTEX);
        int nPoints = displayedPoints.length;
        for (int n = 0; n < nPoints; n++)
        {
            gl.glPushName(n);
            StsGLDraw.drawPoint2d(displayedPoints[n], gl, 6);
            gl.glPopName();
        }
    }


    public void pickMutes(GL gl, StsVelocityProfile velocityProfile)
    {
        gl.glInitNames();
        gl.glPushName(StsVelocityAnalysisEdit.TYPE_GATHER_VERTEX);
        gl.glPushName(StsGather.TOP_MUTE_INDEX);
        StsGLDraw.drawPoint2d(velocityProfile.getTopMute(), gl, 6);
        gl.glPopName();
        gl.glPushName(StsGather.BOT_MUTE_INDEX);
        StsGLDraw.drawPoint2d(velocityProfile.getBottomMute(), gl, 6);
        gl.glPopName();
    }

    /**
     * Note that pushName/popName are outside the glBegin/glEnd loop (required),
     * which makes this more expensive than drawing which is the reason for the
     * separate pick method.
     */
    public void pickLine(GL gl, StsPoint[] semblancePoints)
    {
        gl.glInitNames();
        gl.glPushName(StsVelocityAnalysisEdit.TYPE_GATHER_EDGE);
        int nPoints = semblancePoints.length;
        for (int n = 0; n < nPoints - 1; n++)
        {
            gl.glPushName(n);
            gl.glBegin(GL.GL_LINE_STRIP);
            gl.glVertex2fv(semblancePoints[n].v, 0);
            gl.glVertex2fv(semblancePoints[n + 1].v, 0);
            gl.glEnd();
            gl.glPopName();
        }
    }

    public void restoreVertexPick(StsVelocityProfile velocityProfile)
    {
        if (pickedIndex == -1) return;
        velocityProfile.setProfilePoint(pickedIndex, pickedPoint);

        pickedIndex = -1;
    }

    /**
     * Method for default model editing
     */
    public StsObject defaultEdit(StsMouse mouse)
    {
        try
        {
            int currentButton = mouse.getCurrentButton();
            // middle mouse pick on vertex or profile line results in yes-no dialog asking user if he wants to delete it
            if (currentButton == StsMouse.MIDDLE)
            {
                if (glPanel3d.actionManager.getCurrentAction() instanceof StsVelocityAnalysisEdit)
                {
                    if (currentVelocityProfile == null) return null;

                    int buttonState = mouse.getButtonStateCheckClear(StsMouse.MIDDLE);
                    if (buttonState == StsMouse.PRESSED )
                    {
                        StsMethod pickMethod = new StsMethod(this, "pickOnSeismicView", new Object[]{glPanel, currentVelocityProfile});
                        StsJOGLPick.pick3d(glPanel, pickMethod, StsMethodPick.PICKSIZE_EXTRA_LARGE, StsMethodPick.PICK_ALL);

                        if (StsJOGLPick.hits > 0)
                        {
                            StsPickItem pickItem = StsJOGLPick.pickItems[0];
                            int typePicked = pickItem.names[0];
                            int indexPicked = pickItem.names[1];
                            if (typePicked == StsVelocityAnalysisEdit.TYPE_GATHER_VERTEX)
                            {
                                boolean delete = true; //StsMessage.questionValue(glPanel3d.window, "Delete selected vertex?"); //Processors want it to just delete!!
                                return deletePickedPoint(indexPicked, delete);
                            } else if (typePicked == StsVelocityAnalysisEdit.TYPE_GATHER_EDGE && indexPicked >= 0)
                            {
                                boolean delete = StsYesNoDialog.questionValue(glPanel3d.window, "Delete profile?");
                                return deleteCurrentVelocityProfile(delete);
                            }
                        }
                        mouse.clearButtonState(StsMouse.MIDDLE);
                    }
                }
            } else if (currentButton == StsMouse.LEFT)
            {
                StsPreStackVelocityModel velocityModel = lineSet.getVelocityModel();
                if (velocityModel == null) return null;

                if (currentVelocityProfile == null) return null;
                int currentRow = currentVelocityProfile.row;
                int currentCol = currentVelocityProfile.col;
                StsSemblanceDisplayProperties semblanceDisplayProperties = lineSet.getSemblanceDisplayProperties();
                int ilineDisplayInc = semblanceDisplayProperties.getIlineDisplayInc();
                int xlineDisplayInc = semblanceDisplayProperties.getXlineDisplayInc();
                StsVelocityProfile pickedProfile = StsVelocityProfile.pickNeighborProfiles(currentRow, currentCol, glPanel3d, ilineDisplayInc, xlineDisplayInc, lineSet);
                if (pickedProfile != null)
                {
                    lineSet.jumpToRowCol(new int[]{pickedProfile.row, pickedProfile.col}, glPanel3d.window);
                    return pickedProfile;
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("Exception in StsModel.defaultEdit()", e, StsException.WARNING);
        }
        return null;
    }

    // Change between pick and display modes.
    public void displayMode()
    {
        // ToDo Tom: Put in Display mode on press of D key
    }

    public void editMode()
    {
        // ToDo Tom: Put in Edit mode on press of E key
    }

    public void next()
    {
        lineSet.nextProfile(glPanel3d.window);
        glPanel3d.window.resetFamilyViews();
    }

    public void previous()
    {
        lineSet.previousProfile(glPanel3d.window);
        glPanel3d.window.resetFamilyViews();
    }

    private StsVelocityProfile deletePickedPoint(int indexPicked, boolean delete)
    {
        if (!delete) return currentVelocityProfile;

        if (indexPicked < 0 || indexPicked > currentVelocityProfile.getProfilePoints().length)
            return currentVelocityProfile;

        currentVelocityProfile.deleteProfilePoint(indexPicked);
        superGather.velocityProfileChanged(currentVelocityProfile);

        StsPreStackVelocityModel velocityModel = lineSet.getVelocityModel();
        if (velocityModel != null)
        {
            velocityModel.recomputeCorridor(); // do we need to recompute entire corridor or just a section ?
            if (indexPicked >= 0)
                velocityModel.updateVelocityProfile(currentVelocityProfile, indexPicked);
        }

        viewObjectChangedAndRepaint(this, currentVelocityProfile);
        //        model.viewObjectRepaint(lineSet); // changes 3d view of 3d data (lineSet is cursor3d displayable object)
        //        model.viewObjectRepaint(currentVelocityProfile); // redraws 2d line in 3d view

        //        model.repaintViews(StsView3d.class);
        return currentVelocityProfile;
    }


    private void viewObjectChangedAndRepaint(Object source, StsVelocityProfile velocityProfile)
    {
        /* if(lineSet instanceof StsPreStackLineSet2d)
        {
            StsPreStackLine2d currentLine2d = (StsPreStackLine2d)lineSet.getDataLine(velocityProfile.row, velocityProfile.col);
            model.viewObjectChangedAndRepaint(source, currentLine2d);
        }
        else
        {*/
        model.viewObjectChangedAndRepaint(source, velocityProfile);
        //}
    }

    private StsVelocityProfile deleteCurrentVelocityProfile(boolean delete)
    {
        if (!delete) return currentVelocityProfile;
        currentVelocityProfile.delete();
        superGather.reinitialize();
        viewObjectChangedAndRepaint(this, currentVelocityProfile);
        //        model.viewObjectRepaint(lineSet);
        //        model.viewObjectRepaint(currentVelocityProfile);
        return currentVelocityProfile;
    }

    public void stretch(StsMouse mouse, int index, float factor, byte axis)
    {
        float[][] oldAxisRange = StsMath.copyFloatArray(axisRanges);
        super.stretch(mouse, index, factor);
        recomputeZoomedScaling(oldAxisRange, axis);
    }

    public void stretch(StsMouse mouse, int index, float factor)
    {
        stretch(mouse, index, factor, BOTH);
    }

    /* Calculate the vertical and horizontal scaling after a zoom operation
    */
    public void recomputeZoomedScaling(float[][] oldRanges, byte axis)
    {
        if (lineSet == null) return;
        float inches, axisRange;
        if ((axis == BOTH) || (axis == VERTICAL))
        {
            axisRange = (axisRanges[1][0] - axisRanges[1][1]) / 1000.0f;  //seconds
            inches = getInsetHeight() / (float) getPixelsPerInch();  // inches
            lineSet.getWiggleDisplayProperties().setInchesPerSecond(inches / axisRange);
        }

        if ((axis == BOTH) || (axis == HORIZONTAL))
        {
            float tracesPerInch = lineSet.getWiggleDisplayProperties().getTracesPerInch();
            inches = (float) glPanel3d.getWidth() / (float) getPixelsPerInch();  // inches
            float oldRange = (oldRanges[0][1] - oldRanges[0][0]);
            float unitsPerTrace = (oldRange / inches) / tracesPerInch;
            axisRange = (axisRanges[0][1] - axisRanges[0][0]);  // distance
            float nTraces = axisRange / unitsPerTrace;
            lineSet.getWiggleDisplayProperties().setTracesPerInch(nTraces / inches);
        }
        // Not optimal but without it a zoom operation will not be persisted......SAJ
        // scale parameters now moved to StsClass level and are persisted with propertiesPersistManager. TJL 3/5/07
        //        lineSet.getWiggleDisplayProperties().saveState();
        //        lineSet.getWiggleDisplayProperties().commitChanges();
    }
    /* Calculate the vertical axis limits based on user input
    */

    public void rescaleInchesPerSecond()
    {
        // with our current size, how many samples can we fit ?
        float inches = (float) getInsetHeight() / (float) getPixelsPerInch();
        float nSeconds = inches / lineSet.getWiggleDisplayProperties().getInchesPerSecond();
        //        axisRanges[1][1] = axisRanges[1][1];   do nothing statement
        axisRanges[1][0] = axisRanges[1][1] + (nSeconds * 1000.0f);
        if (axisRanges[1][0] > totalAxisRanges[1][0])
            axisRanges[1][0] = totalAxisRanges[1][0];
        //        glPanel3d.viewChanged = true;
    }

    public void setVertexPick(int pickedIndex, StsVelocityProfile velocityProfile)
    {
        this.pickedIndex = pickedIndex;
        if (pickedIndex == -1)
            pickedPoint = null;
        else
            pickedPoint = new StsPoint(velocityProfile.getProfilePoint(pickedIndex));
    }

    /**
     * clears semblanceView; calling method should execute gl.glEnable(GL.GL_DEPTH_TEST); gl.glEnable(GL.GL_LIGHTING); in a final block
     */
    protected void clearView()
    {
        gl.glDrawBuffer(GL.GL_BACK);
        // glPanel3d.setClearColor(model.project.getBackgroundStsColor());
        // glPanel3d.applyClearColor();
        gl.glDisable(GL.GL_DEPTH_TEST);
        gl.glDisable(GL.GL_LIGHTING);
        clearToBackground(GL.GL_COLOR_BUFFER_BIT);
        // gl.glClear(GL.GL_COLOR_BUFFER_BIT);
    }

    public StsColor getBackgroundColor()
    {
        return lineSet.getWiggleDisplayProperties().getWiggleBackgroundColor();
    }

    protected void displayAxes()
    {
        try
        {
            resetViewPort();
            gl.glDisable(GL.GL_LIGHTING);
            gl.glMatrixMode(GL.GL_PROJECTION);
            gl.glPushMatrix();
            gl.glLoadIdentity();
            glu.gluOrtho2D(0, glPanel3d.getWidth(), 0, glPanel3d.getHeight());
            gl.glMatrixMode(GL.GL_MODELVIEW);

            String axisLabel = getXAxisLabel();
            // Title needs to be Seismic Format not Crossplot
            //				float lineNum = lineSet.getRowNumFromRow(currentRow);
            //				float xlineNum = lineSet.getColNumFromCol(currentCol);
            String titleLabel = null;
            /*                if (is2dSeismic)
            {
                titleLabel = lineSet.getTitleString2d(currentRow, currentCol);
            }
            else*/
            titleLabel = superGather.getGatherDescription();

            //titleLabel = " Inline " + lineNum + " Crossline " + xlineNum;

            drawHorizontalAxis(titleLabel, axisLabel, getScaledHorizontalAxisRange(), model.getProject().getShow2dGrid());
            //                drawHorizontalAxis( gl, titleLabel, axisLabel, axisRanges[0], lineSet.getDisplayGridLines());
            axisLabel = getYAxisLabel();
            //            if(debug) System.out.println("Vertical Axis Range: " + axisRanges[1][0] + " - " + axisRanges[1][1]);
            drawVerticalAxis(axisLabel, getVerticalAxisRange(), model.getProject().getShow2dGrid());

            gl.glMatrixMode(GL.GL_PROJECTION);
            gl.glPopMatrix();
            gl.glMatrixMode(GL.GL_MODELVIEW);

            //            insetViewPort(); // sets up viewport for picking
        }
        catch (Exception e)
        {
            StsException.outputException("StsViewSemblance.display() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public String getXAxisLabel()
    {
        return axisLabels[0];
    }

    public String getYAxisLabel()
    {
        return axisLabels[1];
    }

    /**
     * Override if isVisible axis is scaled from actual axis
     */
    public float[] getScaledHorizontalAxisRange()
    {
        return getHorizontalAxisRange();
    }

    /**
     * Override if isVisible axis is scaled from actual axis
     */
    public float[] getScaledVerticalAxisRange()
    {
        return getVerticalAxisRange();
    }

    public void reshape(GLAutoDrawable drawable, int x, int y, int w, int h)
    {
    }

    protected void displayProfiles()
    {
        StsSemblanceDisplayProperties semblanceDisplayProperties = lineSet.getSemblanceDisplayProperties();
        StsPreStackVelocityModel velocityModel = lineSet.velocityModel;
        int currentRow = superGather.superGatherRow;
        int currentCol = superGather.superGatherCol;
        if (velocityModel != null)
        {
            // display velocity trace from velocity cube if available
            if (semblanceDisplayProperties.showAvgProfile)
            {
                float[][] traceVelocities = velocityModel.getVolumeVelocityProfile(currentRow, currentCol);
                if (traceVelocities != null)
                    StsVelocityProfile.displayTraceVelocity(glPanel3d, traceVelocities, StsColor.BLUE);
                StsPoint[] interpolatedPoints = velocityModel.getInterpolatedVelocityProfilePoints(currentRow, currentCol);
                if (interpolatedPoints != null)
                    velocityModel.drawInterpolatedProfile(glPanel3d, interpolatedPoints, StsColor.GREEN, semblanceDisplayProperties.corridorPercentage);
            }

            if (semblanceDisplayProperties.showLimitProfiles && velocityModel.getNProfiles() > 1)
            {
                //   glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);
                velocityModel.drawProfileCorridor(gl, semblanceDisplayProperties.corridorPercentage);
                //    glPanel3d.resetViewShift(gl);
            }

            // display neighbor profiles if parameters set for them
            if (semblanceDisplayProperties.showNeighborProfiles)
                StsVelocityProfile.displayNeighborProfiles(currentRow, currentCol, glPanel3d, semblanceDisplayProperties.ilineDisplayInc,
                        semblanceDisplayProperties.xlineDisplayInc, semblanceDisplayProperties.showLabels, lineSet, pickedIndex, false);

            currentVelocityProfile = getVelocityProfile(currentRow, currentCol);
            // display picked velocity profile
            if (currentVelocityProfile != null) // && !currentVelocityProfile.isInterpolated())
            {
                // StsVelocityProfile.debugPrintProfile(this, "displayProfiles", currentVelocityProfile);
                if (semblanceDisplayProperties.showIntervalVelocities)  //point labels should plot on top of interval vel lines
                {
                    StsVelocityProfile.displayIntervalVelocities(glPanel3d, superGather, currentVelocityProfile);
                }
                //boolean labelBackground = false;
                //                if(this instanceof StsViewStacks) //black background makes it so you can't see stacks SWC  - 5/12/2009
                //                    labelBackground = true;
                velocityModel.drawProfile(glPanel3d, currentVelocityProfile.getProfilePoints(), StsColor.RED, semblanceDisplayProperties.corridorPercentage, true, true, semblanceDisplayProperties.showLabels, pickedIndex, false);      //reenable wrw
//                StsVelocityProfile.displayOnSemblance(currentVelocityProfile, glPanel3d, StsColor.RED, true, semblanceDisplayProperties.showLabels, pickedIndex, false);  //velocityModel.drawProfile() is unnecessary and not updating correctly with "snapped" picks (pickTime - prevPick < pickThreshold) SWC 8/20/09
                //                StsVelocityProfile.displayProfile(glPanel3d, superGather, currentVelocityProfile);
                if (semblanceDisplayProperties.showIntervalVelocities)
                {
                    StsVelocityProfile.displayIntervalVelocityLables(glPanel3d, superGather, currentVelocityProfile);
                }
            }
        } else
            currentVelocityProfile = null;
    }

    public StsVelocityProfile getVelocityProfile(int row, int col)
    {
        if (row < 0 || col < 0) return null; //causing index out of bounds when loading saved 2d project
        if (lineSet.velocityModel == null) return null;
        return lineSet.velocityModel.getExistingVelocityProfile(row, col);
    }

    /**
     * Display the current crossplot view
     *
     * @param component the drawable component
     */
    public void display(GLAutoDrawable component)
    {
        if (glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }
        if (debug) debugThread("entered display()");

        if (lineSet == null || axisRanges == null)
        {
            clearView();
            return;
        }

        boolean displayAxis = lineSetClass.getDisplayAxis();
        setInsets(displayAxis);
        if (superGather == null) superGather = lineSet.getSuperGather(glPanel3d.window);
        if (isDataUnavailable())
        {
            checkViewChanged();
            if (debug) debugThread("clearing screen: data not available");
            clearView();
            if (displayAxis) displayAxes();
            return;
        }
        if (checkComputeData())
        {
            // if(debug) debugThread("clearing screen: computing data");
            // clearView();
            // if(displayAxis) displayAxes();
            return;
        }

        if (isDataChanged())
        {
            if (isDataNull())
            {
                if (debug) debugThread("clearing screen: isDataChangedAndNull");
                checkViewChanged();
                clearView();
                displayAxes();
                return;
            }
        }

        checkTextureSizeAndRange();

        if (volumeChanged)
        {
            computeProjectionMatrix();
            computeModelViewMatrix();
            volumeChanged = false;
        }

        if (glPanel3d.viewPortChanged)
        {
            computePixelScaling();
            glPanel3d.resetViewPort();
            glPanel3d.viewPortChanged = false;
        }

        if (displayAxis) insetViewPort();
        checkViewChanged();

        try
        {
            if (cursorButtonState != StsMouse.CLEARED)
            {
                // if this is window where cursor is being dragged and we have focus, draw foreground cursor.
                // If not the window where cursor is being dragged, but we are displaying cursor here,
                // draw the windows;
                if (isCursorWindow && glPanel3d.hasFocus() || !isCursorWindow)
                {
                    //                    System.out.println("draw cursor on Semblance");
                    drawForeground(gl);
                    if (cursorButtonState != StsMouse.CLEARED) return;
                }
                if (!isCursorWindow && glPanel3d.hasFocus())
                {
                    ;
                }
            }

            clearView();
            displayData(gl, glu);

            if (displayAxis)
            {
                displayAxes();
                insetViewPort(); // sets up viewport for picking
            }
            displayProfiles();
            StsVelocityCursor.getStsVelocityCursor().displayCursor(gl, this);
        }
        catch (Exception e)
        {
            exceptionCleanup(glPanel3d);
            StsException.outputWarningException(this, "display", e);
        }
        finally
        {
            gl.glEnable(GL.GL_DEPTH_TEST);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    private void checkViewChanged()
    {
        if (!glPanel3d.viewChanged) return;
        computeProjectionMatrix();
        computeModelViewMatrix();
        glPanel3d.viewChanged = false;
    }

    protected void initializeView()
    {
        super.initializeView();
        matchLockedWindows();
    }

    // methods used by display() for all prestackSeismic displays; override as needed in concrete subclasses
    protected void debugThread(String message)
    {
        System.out.println(Thread.currentThread().getName() + " window " + glPanel3d.window.getTitle() + " view class " + StsToolkit.getSimpleClassname(this) + " " + message);
    }

    protected boolean checkComputeData()
    {
        return false;
    }

    protected boolean isDataChanged()
    {
        return false;
    }

    protected boolean isDataNull()
    {
        return false;
    }

    protected boolean isDataUnavailable()
    {
        this.getWindow();
        // boolean isLeftMouseDown = isCursor3dDragging();
        int currentRow = superGather.superGatherRow;
        int currentCol = superGather.superGatherCol;
        return superGather == null || !superGather.filesOk || currentRow == -1 || currentCol == -1;
    }

    protected void checkTextureSizeAndRange()
    {
    }

    protected void initializeAxisLabels()
    {
    }

    protected void displayData(GL gl, GLU glu)
    {
    }

    protected void exceptionCleanup(StsGLPanel3d glPanel3d)
    {
    }

    public void drawForeground(GL gl)
    {
    }

    public void resetToOrigin()
    {
        float xRange = axisRanges[0][1] - axisRanges[0][0];
        float yRange = axisRanges[1][0] - axisRanges[1][1];
        axisRanges[0][0] = totalAxisRanges[0][0];
        axisRanges[0][1] = totalAxisRanges[0][0] + xRange;
        axisRanges[1][0] = totalAxisRanges[1][1] + yRange;
        axisRanges[1][1] = totalAxisRanges[1][1];
        computePixelScaling();
        repaint();
    }

    public Class getDisplayableClass()
    {
        return StsPreStackLineSet.class;
    }

    public String getName()
    {
        return getShortViewName() + " of " + lineSet.getName();
    }

    /**
     * sets currentVelocityProfile
     *
     * @param profile
     */
    public void setVelocityProfile(StsVelocityProfile profile)
    {
        currentVelocityProfile = profile;
    }

    /**
     * returns currentVelocityProfile
     *
     * @return
     */
    public StsVelocityProfile getVelocityProfile()
    {
        return currentVelocityProfile;
    }
}