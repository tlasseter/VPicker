package com.Sts.MVC;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.util.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author T.J.Lasseter
 * @version c63d
 */

/** Class defining a prefix or family of windows. Each family will have one parent and any number of children.
 * The parent will control all the children and is equipped with a 3D Cursor Panel. The default window (StsWin3d)
 * is a parent of the default window family.
 */

public class StsWindowFamily implements /*Serializable,*/ StsSerializable
{
	//StsWin3dBase parent;
	StsWin3dBase[] windows;
    StsViewTimeSeries[] timeWindows;
    StsViewGauges[] gaugeWindows;

	/**
	 * StsWindowFamily constructor for default window (StsWin3d). This parent has all decorations including 3D cursor
	 * panel, status panel and object and workflow tree.
	 **/
	public StsWindowFamily()
	{
		// null const
    }

	 public StsWindowFamily (StsWin3d parent)
	{
		//this.parent = parent;
		windows = new StsWin3dBase[] { parent };
	}

	/**
	 * StsWindowFamily constructor for full function window (StsWin3dFull). This parent does not include the
	 * object and workflow panel.
	 * @param parent the StsWin3dFull
	 */
	public StsWindowFamily(StsWin3dFull parent)
	{
		//this.parent = parent;
		windows = new StsWin3dBase[] { parent };
	}

	/**
	 * Add a window to this family
	 * @param window the auxillary window to add
	 */
	public void addWindow(StsWin3dBase window)
	{
		windows = (StsWin3dBase[])StsMath.arrayAddElement(windows, window);
	}
    public void addTimeSeriesWindow(StsViewTimeSeries window)
    {
        timeWindows = (StsViewTimeSeries[])StsMath.arrayAddElement(timeWindows, window);
    }
    public void addGaugeWindow(StsViewGauges window)
    {
        gaugeWindows = (StsViewGauges[])StsMath.arrayAddElement(gaugeWindows, window);
    }
	/**
	 * Delete a window from this family
	 * @param auxWindow the auxillary window to delete
	 */
	public void deleteAuxWindow(StsWin3dBase auxWindow, int familyIndex)
	{
		windows = (StsWin3dBase[])StsMath.arrayDeleteElement(windows, auxWindow);
		renumberWindows(familyIndex);
	}

    public void deleteAuxTimeWindow(StsViewTimeSeries auxWindow, int familyIndex)
    {
        timeWindows = (StsViewTimeSeries[])StsMath.arrayDeleteElement(timeWindows, auxWindow);
        renumberTimeWindows(familyIndex);
	}

    public void deleteAuxGaugeWindow(StsViewGauges auxWindow, int familyIndex)
    {
        gaugeWindows = (StsViewGauges[])StsMath.arrayDeleteElement(gaugeWindows, auxWindow);
        renumberGaugeWindows(familyIndex);
	}

	public void renumberWindows(int familyIndex)
	{
		String familyName;
		if(familyIndex == 0)
			familyName = "Main Group";
		else
		    familyName = "Group " + familyIndex;
		for(int w = 1; w < windows.length; w++)
		{
			StsWin3dBase auxWindow = (StsWin3dBase)windows[w];
			windows[w].setTitle( new String("   Auxiliary Window, " + familyName + ", Window " + w));
		}
	}

    public void renumberTimeWindows(int familyIndex)
    {
        String familyName;
        if(familyIndex == 0)
            familyName = "Main Group";
        else
            familyName = "Group " + familyIndex;
        for(int w = 1; w < timeWindows.length; w++)
        {
            StsViewTimeSeries auxWindow = (StsViewTimeSeries)timeWindows[w];
            timeWindows[w].setTitle( new String("   Auxiliary Time Window, " + familyName + ", Window " + w));
        }
    }


    public void renumberGaugeWindows(int familyIndex)
    {
        String familyName;
        if(familyIndex == 0)
            familyName = "Main Group";
        else
            familyName = "Group " + familyIndex;
        for(int w = 1; w < gaugeWindows.length; w++)
        {
            StsViewGauges auxWindow = (StsViewGauges)gaugeWindows[w];
            gaugeWindows[w].setTitle( new String("   Auxiliary Gauge Window, " + familyName + ", Window " + w));
        }
    }
    /**
	 * Get the number of windows in this family, including the parent
	 * @return total of windows in family
	 */
	public int getNumberWindows()
	{
		return windows.length;
	}
    public int getNumberTimeWindows()
    {
        return timeWindows.length;
	}
    public int getNumberGaugeWindows()
    {
        return gaugeWindows.length;
	}
	/**
	 * Set all the windows in this family to the default view
	 */
    /*
    public void setDefaultView()
	{
		for(int n = 0; n < windows.length; n++)
			windows[n].getSingleViewGlPanel3d().setDefaultView();
	}
    */
	/**
	 * Redisplay all windows in this family
	 */
	public void win3dDisplay()
	{
        Iterator<StsView> windowViewIterator = getWindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = windowViewIterator.next();

			// jbw restore problem -- panel may not exist yet
			if (view.glPanel3d != null)
				view.glPanel3d.repaint();
		}
        if(timeWindows != null)
        {
            for(int n = 0; n < timeWindows.length; n++)
            {
                timeWindows[n].repaint();
            }
        }
        if(gaugeWindows == null) return;
        for(int n = 0; n < gaugeWindows.length; n++)
        {
            gaugeWindows[n].repaint();
        }
	}

	public void viewObjectRepaint(Object source, Object object)
	{
        WindowViewIterator windowViewIterator = new WindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = (StsView)windowViewIterator.next();
            view.viewObjectRepaint(source, object);
		}

        // Time Series Views
        if(timeWindows != null)
        {
            for(int n = 0; n < timeWindows.length; n++)
                timeWindows[n].viewObjectRepaint(source, object);
        }

        // Gauge Views
        if(gaugeWindows == null) return;
        for(int n = 0; n < gaugeWindows.length; n++)
            gaugeWindows[n].viewObjectRepaint(source, object);
	}

	public void repaintViews(Class viewClass)
	{
        WindowViewIterator windowViewIterator = new WindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = windowViewIterator.next();
                if(viewClass.isAssignableFrom(view.getClass()))
                    view.glPanel3d.repaint();
        }
    }

	public void repaintViews()
	{
        WindowViewIterator windowViewIterator = new WindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = (StsView)windowViewIterator.next();
            view.glPanel3d.repaint();
        }
    }

    public void viewChanged(Class viewClass)
	{
        WindowViewIterator windowViewIterator = new WindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = (StsView)windowViewIterator.next();
                if(viewClass.isAssignableFrom(view.getClass()))
                    view.viewChanged();
        }
    }

    public void viewChanged()
	{
        WindowViewIterator windowViewIterator = new WindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = (StsView)windowViewIterator.next();
            view.viewChanged();
        }
    }
    /** Notify views that properties of an object that they may be viewing is changed.
	 *  If they are interested in this object, they return true in which case repaint is called on glPanel3d.
	 */
	public boolean viewObjectChanged(Object source, Object object)
	{
        boolean changed = false;
        WindowViewIterator windowViewIterator = new WindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = windowViewIterator.next();
			changed = changed | view.viewObjectChanged(source, object);
		}
        // Time Series Views
        if(timeWindows != null)
        {
            for(int n = 0; n < timeWindows.length; n++)
                changed = changed | timeWindows[n].viewObjectChanged(source, object);
        }
        // Gauge Views
        if(gaugeWindows == null) return changed;
        for(int n = 0; n < gaugeWindows.length; n++)
            changed = changed | gaugeWindows[n].viewObjectChanged(source, object);
        return changed;
    }

    public boolean viewObjectChangedAndRepaint(Object source, Object object)
    {
        boolean changed = false;
        WindowViewIterator windowViewIterator = new WindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = (StsView)windowViewIterator.next();
            if(view.viewObjectChanged(source, object))
            {
                changed = true;
                viewObjectRepaint(source, object);
            }
        }
        // Time Series Views
        if(timeWindows != null)
        {
            for(int n = 0; n < timeWindows.length; n++)
            {
                if(timeWindows[n].viewObjectChanged(source, object))
                {
                     changed = true;
                    timeWindows[n].viewObjectRepaint(source, object);
                }
            }
        }
        // Gauge Views
        if(gaugeWindows == null) return changed;
        for(int n = 0; n < gaugeWindows.length; n++)
        {
            if(gaugeWindows[n].viewObjectChanged(source, object))
            {
                 changed = true;
                 gaugeWindows[n].viewObjectRepaint(source, object);
            }
        }
        return changed;
    }
    /**
	 * Set the background color for all windows in this family
	 * @param color the color to change to
	 */
/*
    public void setClearColor(StsColor color)
	{
        WindowViewIterator windowViewIterator = new WindowViewIterator();
        while(windowViewIterator.hasNext())
        {
            StsView view = (StsView)windowViewIterator.next();
            StsGLPanel3d glPanel3d = view.glPanel3d;
            glPanel3d.setClearColor(color);
            glPanel3d.repaint();
        }
    }
*/
    public Iterator<StsWin3dBase> getWindowIterator()
    {
        return new WindowIterator();
    }

    private class WindowIterator implements Iterator<StsWin3dBase>
    {
        int nWindows;
        int nWindow = 0;
        StsWin3dBase window;

        public WindowIterator()
        {

            nWindows = (windows == null ? 0 : windows.length);
        }

        public boolean hasNext()
        {
            if(nWindow >= nWindows) return false;
            window = windows[nWindow++];
            return true;
        }

        public StsWin3dBase next()
        {
            return window;
        }

        public void remove()
        {
        }
    }

    public Iterator<StsView> getWindowViewIterator()
    {
        return new WindowViewIterator();

    }

    private class WindowViewIterator implements Iterator<StsView>
    {
        Iterator<StsWin3dBase> windowIterator;
        Iterator<StsView> viewIterator = null;

        public WindowViewIterator()
        {
        }

        public boolean hasNext()
        {
            if(windowIterator == null)
            {
                windowIterator = getWindowIterator();
                if(windowIterator.hasNext())
                {
                    StsWin3dBase window = windowIterator.next();
                    viewIterator = window.getViewIterator();
                }
            }
            else
            {
                if(viewIterator.hasNext()) return true;
                if(windowIterator.hasNext())
                {
                    StsWin3dBase window = windowIterator.next();
                    viewIterator = window.getViewIterator();
                }
            }
            if(viewIterator == null) return false;
            return viewIterator.hasNext();
        }

        public StsView next()
        {
            if(viewIterator == null) return null;
            return viewIterator.next();
        }

        public void remove()
        {
        }
    }

    public Iterator<StsView> getWindowViewIteratorOfType(Class viewType)
    {
        return new WindowViewIteratorOfType(viewType);
    }

   private class WindowViewIteratorOfType extends WindowViewIterator
    {
        Class viewType;

        public WindowViewIteratorOfType(Class viewType)
        {
            this.viewType = viewType;
        }

        public boolean hasNext()
        {
            if(windowIterator == null)
            {
                windowIterator = getWindowIterator();
                if(windowIterator.hasNext())
                {
                    StsWin3dBase window = windowIterator.next();
                    viewIterator = window.getViewIteratorOfType(viewType);
                }
            }
            else
            {
                if(viewIterator.hasNext()) return true;
                if(windowIterator.hasNext())
                {
                    StsWin3dBase window = windowIterator.next();
                    viewIterator = window.getViewIteratorOfType(viewType);
                }
            }
            return viewIterator.hasNext();
        }
    }

    /**
	 * Adjust the cursor for all windows in this family to new direction and position
	 * @param dirCoor new position
	 */
	public void adjustCursor(int dir, float dirCoor)
	{
        for(int n = 0; n < windows.length; n++)
			windows[n].adjustCursor(dir, dirCoor);
	}

    public void adjustCursorXY(float xCoor, float yCoor)
    {
        for(int n = 0; n < windows.length; n++)
            windows[n].adjustCursorXY(xCoor, yCoor);
    }

	public void adjustCursorAndSlider(int dir, float dirCoor)
	{
        getParent().adjustCursorAndSlider(dir, dirCoor);
	}

	/**
	 * Enable the cursor in the specified window
	 * @param currentDirection the direction to enable
	 * @param enable enable or disable
	 */
	public void setSelectedDirection(int currentDirection, boolean enable)
	{
        StsWin3dFull win3d = (StsWin3dFull)windows[0];
        win3d.getCursor3d().setSelectedDirection(currentDirection, enable);
	}

	public void incrementCursor()
	{
        StsWin3dFull win3d = (StsWin3dFull)windows[0];
        win3d.incrementCursor(win3d.getCursor3d().getCurrentDirNo());
        win3d.model.win3dDisplayAll();        
	}

    public void decrementCursor()
	{
        StsWin3dFull win3d = (StsWin3dFull)windows[0];
        win3d.decrementCursor(win3d.getCursor3d().getCurrentDirNo());
        win3d.model.win3dDisplayAll();
	}
	/**
	 * Set the current direction of the cursor in the current window
	 * @param currentWindow the specified window
	 * @param currentDirection the direction to enable
	 */
	public void setCursor3dCurrentDirNo(StsWin3dBase currentWindow, int currentDirection)
	{
		for(int n = 0; n < windows.length; n++)
		{
			windows[n].getCursor3d().setCurrentDirNo(currentDirection);
		}
	}

    public void reinitializeCursor3d()
    {
        Iterator windowIterator = getWindowIterator();
        while(windowIterator.hasNext())
        {
            StsWin3dBase window = (StsWin3dBase)windowIterator.next();
            window.getCursor3d().initialize();
            window.getCursor3d().resetInitialCursorPositions();
            window.setDefaultView();
        }
    }

    /**
	 * Clear the cursor data (subImage2D) display for all windows in the family.
	 * Will require that the subImage2D with new data be reconstructed on the next display.
	 */
    /*
    public void clearDataDisplays()
	{
		for(int n = 0; n < windows.length; n++)
		{
			StsWin3dBase auxWindow = windows[n];
			if(auxWindow.getGlPanel3d() == null) return;
			StsView[] views = auxWindow.getGlPanel3d().views;
			if(views == null) continue;
			for(int v = 0; v < views.length; v++)
			{
				views[v].clearDataDisplay();
			}
		}
	}
    */

    /**
	 * Clear the cursor texture for all windows.
	 * Will require that the background texture and subImage2D be reconstructed on the next display.
	 */
    /*
    public void clearTextureDisplays()
	{
		for(int w = 0; w < windows.length; w++)
		{
            StsWin3dBase window = windows[w];
			if(window.getGlPanel3d() == null) continue;
            window.getGlPanel3d().clearAllTextures();
		}
	}
    */

    public void deleteClassTextureDisplays(Class displayableClass)
	{
        StsCursor3d cursor3d = getCursor3d();
        if(cursor3d == null) return;
        cursor3d.deleteClassTextureDisplays(displayableClass);
        // StsTextureList.deleteViewClassTextures(viewClass);
    }

    public void deleteCursor3dTextures(StsObject object)
	{
        for(StsWin3dBase window : windows)
        {
            StsCursor3d cursor3d = window.getCursor3d();
            if(cursor3d != null) cursor3d.deleteCursor3dTextures(object);
		}
	}

	/**
	 * Crop has changed: adjust texture coordinates on each window.
	 * Will require that the background texture and subImage2D be reconstructed on the next display.
	 */
	public void cropChanged()
	{
		for(int w = 0; w < windows.length; w++)
		{
			StsWin3dBase auxWindow = windows[w];
			auxWindow.getCursor3d().cropChanged();
		}
	}

	public void subVolumeChanged()
	{
        StsCursor3d cursor3d = getCursor3d();
        if(cursor3d == null) return;
        cursor3d.subVolumeChanged();
    }

	/**
	 * Close the windows in this family
	 */
	public void close(GL gl)
	{
        if(gl != null)
        {
            StsCursor3d cursor3d = getCursor3d();
            if(cursor3d != null)
                cursor3d.deleteAllTextures(gl);
        }
        if(windows != null)
        {
            for(int n = windows.length - 1; n >= 0; n--)
            {
                StsWin3dBase window = windows[n];
                window.dispose();
            }
            windows = null;
        }

        // Time series windows
        if(timeWindows != null)
        {
            for(int n = timeWindows.length - 1; n >= 0; n--)
            {
                StsViewTimeSeries window = timeWindows[n];
                window.dispose();
            }
            timeWindows = null;
        }

        // Gauge windows
        if(gaugeWindows != null)
        {
            for(int n = gaugeWindows.length - 1; n >= 0; n--)
            {
                StsViewGauges window = gaugeWindows[n];
                window.dispose();
            }
            gaugeWindows = null;
        }
	}
	/**
	 * Get the parent of this family
	 * @return parent window (StsWin3d or StsWin3dFull)
	 */
	public StsWin3dFull getParent()
	{
		if(windows == null || windows.length == 0) return null;
		return (StsWin3dFull)windows[0];
	}

    public StsCursor3d getCursor3d()
    {
        StsWin3dFull parentWindow = getParent();
        if(parentWindow == null) return null;
        return parentWindow.cursor3d;
    }

    /**
	 * Get the windows in this family including the parent window
	 * @return ArrayList of windows (StsWin3dBase)
	 */
	public StsWin3dBase[] getWindows()
	{
		return windows;
	}

    public StsViewTimeSeries[] getTimeWindows()
    {
        return timeWindows;
	}

    public StsViewGauges[] getGaugeWindows()
    {
        return gaugeWindows;
	}
	public void win3dDeleteAllDisplayLists()
	{
		for(int w = 0; w < windows.length; w++)
		{
			StsWin3dBase auxWindow = windows[w];
			auxWindow.getCursor3d().cropChanged();
		}
	}

	}