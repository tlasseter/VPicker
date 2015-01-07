package com.Sts.UI.DataTransfer;

/**
 * <p>Title: S2S development</p>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p>Copyright: Copyright (c) 2005</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version c62e
 */

/** Interface implemented by objects which are listeners on the StsObjectTransferPanel */
public interface StsObjectTransferFace
{
    /** objects available to be put in left-hand list */
    public Object[] getAvailableObjects();
    /** objects added to the right-hand side */
    public void addObjects(Object[] objects);
    /** objects removed from the right and put back on the left */
    public void removeObjects(Object[] objects);
    /** all objects removed from the right */
    public void removeAllObjects();
    /** object has been single-selected on the right (post info about it for example */
    public void objectSelected(Object selectedObject);
}