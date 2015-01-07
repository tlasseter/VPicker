package com.Sts.MVC.ViewWell;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;
import com.Sts.Actions.*;
import com.magician.fonts.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.GLU;
import java.awt.event.*;

/**
 * <p>Title: S2S development</p>
 * <p/>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p/>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public abstract class StsWellView extends StsView2d
{
    public StsObject displayedObject;
    protected StsWell well;
    transient protected StsWellViewModel wellViewModel;
    /** min mdepth for this view */
	protected double windowMdepthMin;
	/** max mdepth for this view */
	protected double windowMdepthMax;
    transient public StsJPanel curveNameBackPanel;
    /* inner panel for splitpane */
    protected transient StsJPanel innerPanel;

    abstract public boolean initializeView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow);

    abstract public void savePixels(boolean savem);

    public StsWellView()
    {
    }

    public void init(GLAutoDrawable drawable)
    {
        if (isViewGLInitialized) return;
        super.init(drawable);
        glPanel.panelViewChanged = true;
        isViewGLInitialized = true;
    }

    public StsWellView(StsObject displayedObject)
    {
        this.displayedObject = displayedObject;
    }

    public void clearToBackground(int clearBits)
    {
        super.clearToBackground(clearBits);
        StsColor backgroundColor = getBackgroundColor();
        innerPanel.setBackground(backgroundColor.getColor());
    }


    public void viewPortChanged(int x, int y, int width, int height)
    {
        glPanel.viewChanged = true;
    }

    public void computeProjectionMatrix()
    {
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        glu.gluOrtho2D(0, getWidth(), 0, getHeight());
        gl.glMatrixMode(GL.GL_MODELVIEW);
        gl.glLoadIdentity();
    }


    public void computePixelScaling()
    {
    }

    public void moveWindow(StsMouse mouse)
    {
        if (mouse.getButtonStateCheckClear(StsMouse.VIEW) == StsMouse.DRAGGED)
            wellViewModel.moveWindow(mouse);
    }

    public void keyReleased(StsMouse mouse, KeyEvent e)
    {
		int keyCode = KeyEvent.VK_UNDEFINED;
   		if(this.keyCode == keyCode) return;
        if(mouseKeyDebug) StsException.systemDebug(this, "keyReleased", "keyCode changed from " + getKeyText(this.keyCode) + " to " + getKeyText(keyCode));
		this.keyCode = keyCode;
        wellViewModel.keyReleased(this, mouse, e);
    }

	public boolean keyPressed(StsMouse mouse, KeyEvent e)
	{
		int keyCode = e.getKeyCode();
		if(this.keyCode == keyCode) return false;
		if(mouseKeyDebug) StsException.systemDebug(this, "keyPressed", "keyCode changed from " + getKeyText(this.keyCode) + " to " + getKeyText(keyCode));
		this.keyCode = keyCode;
		return wellViewModel.keyPressed(this, mouse, e);
    }
    /*
        public void init(GLAutoDrawable component)
        {
            if(isGLInitialized)return;
            gl.glShadeModel(GL.GL_FLAT);
            gl.glEnable(GL.GL_ALPHA_TEST);
            gl.glAlphaFunc(GL.GL_GREATER, 0.5f);
            computeProjectionMatrix();
            computeModelViewMatrix();
            initializeFont(gl);
            isGLInitialized = true;
        }
    */
    public void performMouseAction(StsActionManager actionManager, StsMouse mouse)
    {
        try
        {
            int but = mouse.getCurrentButton();
            {
                if (but == StsMouse.LEFT)
                {
                    int buttonState = mouse.getButtonState(StsMouse.LEFT);
                    if (buttonState == StsMouse.PRESSED)
                    {
                        wellViewModel.setCursorPosition(mouse.getMousePoint().y, this);
                    }
                    else if (wellViewModel.cursorPicked && buttonState == StsMouse.DRAGGED)
                    {
                        wellViewModel.moveCursor(mouse.getMousePoint().y, this);
                    }
                    else if (buttonState == StsMouse.RELEASED)
                    {
                        StsAction currentAction = actionManager.getCurrentAction();
                        if(currentAction != null) currentAction.performMouseAction(mouse, glPanel);
                        wellViewModel.cursorPicked = false;
                        wellViewModel.savePixels(false);
                        wellViewModel.display();
                    }

                }

                /** If any right mouse action, move view */
                if (but == StsMouse.VIEW && (StsMouse.VIEW != StsMouse.LEFT))
                {
                    wellViewModel.moveWindow(mouse);
                }
// jbw in a remmaped mouse world, mb2 drags the wells
                if (but == StsMouse.MIDDLE && (StsMouse.VIEW == StsMouse.LEFT))
				{
					wellViewModel.moveWindow(mouse);
                }
                /** If middle mouse button clicked, terminate any active function. */

                /** If none active, trigger pop-up menu */
                if (but == StsMouse.POPUP)
                {
                    if(mouse.getButtonStateCheckClear(StsMouse.POPUP) == StsMouse.PRESSED)
                    {
                        showPopupMenu(mouse);
                        //mouse.clearButtonState(StsMouse.MIDDLE, StsMouse.CLEARED);
                    }
					if(mouse.getButtonStateCheckClear(StsMouse.POPUP) == StsMouse.RELEASED)
					{
						cancelPopupMenu(mouse);
						//mouse.clearButtonState(StsMouse.MIDDLE, StsMouse.CLEARED);
					}


                }
            }
        }
        catch (Exception ex)
        {
            StsException.outputWarningException(this, "performMouseAction", ex);
        }
    }

	public void drawCursor(GL gl, int left, int right, int height, boolean displayZ)
	{

		wellViewModel.push2DOrtho(gl, glu, right, height);

		gl.glLineWidth(1.0f);
        getForegroundColor().setGLColor(gl);
		gl.glBegin(GL.GL_LINES);
		gl.glVertex2i(left, wellViewModel.cursorY);
		gl.glVertex2i(right, wellViewModel.cursorY);
		gl.glEnd();

		if(wellViewModel.cursorPicked && displayZ)
		{
			displayZ(gl, left, right);
		}

	    wellViewModel.pop2DOrtho(gl);

	}

    private void displayZ(GL gl, int left, int right)
    {
        String depthLabel = wellViewModel.getMdepthStringFromGLCursorY();
        int locationX = getDisplayZHorizPosition(left, right, depthLabel, horizontalFont);
        int locationY = StsMath.minMax(wellViewModel.cursorY, 0, wellViewModel.displayHeight - wellViewModel.fontHeight);

        //String depthLabel = labelFormatter.format(cursorZ);
        getForegroundColor().setGLColor(gl);
        // StsColor.BLACK.setGLColor(gl);
        StsGLDraw.fontOutput(gl, locationX, locationY, depthLabel, horizontalFont);
//        StsGLDraw.highlightedFontHelvetica12(gl, locationX, locationY, depthLabel);
    }

    protected int getDisplayZHorizPosition(int left, int right, String zLabel, GLBitmapFont font)
    {
        return right;
    }

    public void drawCursorAtMdepth(GL gl, float left, float right, boolean displayZ, StsColor color)
    {
        gl.glLineWidth(1.0f);
        color.setGLColor(gl);
        // float z = (float) getMdepthFromMouseY(displayHeight - cursorY);
        gl.glBegin(GL.GL_LINES);
        float mdepth = (float) wellViewModel.cursorMdepth;
        gl.glVertex2f(left, mdepth);
        gl.glVertex2f(right, mdepth);
        gl.glEnd();

        if (wellViewModel.cursorPicked && displayZ)
        {
            displayZ(gl, left, mdepth);
        }
    }

    private void displayZ(GL gl, float right, float z)
	{
		String depthLabel = wellViewModel.convertMDepthToZString(z);
        getForegroundColor().setGLColor(gl);
        // StsColor.BLACK.setGLColor(gl);
		StsGLDraw.fontHelvetica12(gl, right, z, depthLabel);
	}

	public void drawMarkers(GL gl, GLU glu, int markerLeftY)
	{

		StsObjectRefList markers = well.getMarkers();
		if(markers == null)return;

		int nMarkers = markers.getSize();
		gl.glLineWidth(1.0f);
		int markerLeft = markerLeftY + 2;
		int markerRight = markerLeft + 10;
        int markerCenter = markerLeft + 6;
        double unitsPerPixel = wellViewModel.zoomLevel.unitsPerPixel;
		try
		{
			wellViewModel.initializeFont(gl); // jbw
			for(int n = 0; n < nMarkers; n++)
			{
				StsWellMarker marker = (StsWellMarker)markers.getElement(n);
				float z = marker.getLocation().getM();
				int y = getHeight() - (int)((z - windowMdepthMin) / unitsPerPixel);

				if((marker.getMarker().getType() == StsMarker.SURFACE) || (marker.getMarker().getType() == StsMarker.OFFSET_SURFACE))
				{
					gl.glLineWidth(3.0f);
					marker.getStsColor().setGLColor(gl);
					gl.glBegin(GL.GL_LINES);
					gl.glVertex2i(markerLeft, y);
					gl.glVertex2i(markerRight, y);
					gl.glEnd();
					gl.glLineWidth(1.0f);
				}
				else if(marker instanceof StsPerforationMarker)
				{
					float len = ((StsPerforationMarker)marker).getLength();
					int topY = getHeight() - (int)(((z+(len/2.0f)) - windowMdepthMin) / unitsPerPixel);
					int btmY = getHeight() - (int)(((z-(len/2.0f)) - windowMdepthMin) / unitsPerPixel);
					int midX = markerRight - 5;
					StsColor aColor = new StsColor(marker.getStsColor(), 0.5f);
					aColor.setGLColor(gl);
					gl.glLineWidth(10.0f);
					gl.glEnable(GL.GL_LINE_STIPPLE);
					gl.glLineStipple(1, StsGraphicParameters.dottedLine);
					gl.glBegin(GL.GL_LINE_STRIP);
					gl.glVertex2f(midX, topY);
					gl.glVertex2f(midX, btmY);
					gl.glEnd();
					gl.glDisable(GL.GL_LINE_STIPPLE);
					gl.glLineWidth(1.0f);
				}
				else if(marker instanceof StsEquipmentMarker)
				{
					byte subType = ((StsEquipmentMarker)marker).getSubType();
					StsColor aColor = marker.getStsColor();
					StsGLDraw.drawEquipment2d(markerCenter, y, gl, aColor, subType);
				}
				else if(marker instanceof StsFMIMarker)
				{
					gl.glLineWidth(1.0f);
					marker.getStsColor().setGLColor(gl);
					StsGLDraw.drawCircle2d(new float[] {markerCenter, y}, marker.getStsColor(), gl, StsBitmap.SMALL);
					gl.glBegin(GL.GL_LINES);
					gl.glVertex2i(markerCenter, y-2);
					// Compute tadpole endpoint
					float dX = (float)Math.sin(((StsFMIMarker)marker).getDip()) * 10.0f;
					float dY = (float)Math.cos(((StsFMIMarker)marker).getDip()) * 10.0f;
					gl.glVertex2i((int)(markerCenter - dX), (int)(y - dY));
					gl.glEnd();
					gl.glLineWidth(1.0f);
				}
				else
				{
					marker.getStsColor().setGLColor(gl);
					gl.glBegin(GL.GL_LINES);
					gl.glVertex2i(markerLeft, y);
					gl.glVertex2i(markerRight, y);
					gl.glEnd();
				}

				marker.getStsColor().setGLColor(gl);
				gl.glRasterPos2i(markerRight + 2, y - wellViewModel.fontHeight / 2);
				wellViewModel.font.drawString(marker.getName());
			}
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawMarkers", e);
        }
    }

    public boolean viewObjectRepaint(Object source, Object object)
    {
        return false;
    }
    public boolean viewObjectChanged(Object source, Object object)
    {
        return false;
    }; // { return false; }

}
