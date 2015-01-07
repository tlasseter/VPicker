package com.Sts.Interfaces;

import com.Sts.IO.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 2, 2008
 * Time: 9:43:54 AM
 * To change this template use File | Settings | File Templates.
 */
public interface StsBasicFileTransferFaceNew
{
    void addFiles(StsAbstractFile[] files);
    void removeFiles(StsAbstractFile[] files);
    void removeAllFiles();
    void fileSelected(StsAbstractFile selectedFile);
    void availableFileSelected(StsAbstractFile availableFile);
    boolean hasDirectorySelection();
}