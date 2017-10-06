import requests
from dicom.examples import anonymize
import os

def search_vna(user, pw, accNum=None, study=None, series=None, region='test', limit=None, modality="MR"):
    if region == 'test':
        host = 'vnatest1vt'
        port = '8083'
    elif region == 'prod':
        host = '10.47.11.220'
        port = '8083'
    else:
        raise ValueError("Unsupported region")


    url = ''.join(['http://', host, ':', port,
                   "/AcuoREST/dicomrs/search/studies"])

    if accNum is not None:
        url += "?AccessionNumber=" + accNum

    elif study is not None:
        url += "/" + study + "/series"

        if series is not None:
            url += "/" + series + "/instances"
        elif modality is not None:
            url += "?Modality=" + modality

    if limit is not None:
        if "?" in url:
            url += "&"
        else:
            url += "?"
        url += "limit=" + str(limit)
    #url += "&includefield=all"

    r = requests.get(url, auth=(user, pw)) #, headers=headers
    #if r.status_code != 200:
        #raise ValueError("Invalid request (response code %d) for URL: %s" % (r.status_code, url))
        
    return r, url


def retrieve_vna(user, pw, filename, study=None, series=None, instance=None, region='test', metadata=False):
    """If metadata is true, filename should end in xml. Else end in dcm."""

    if region == 'test':
        host = 'vnatest1vt'
        port = '8083'
    elif region == 'prod':
        host = '10.47.11.220'
        port = '8083'
    else:
        raise ValueError("Unsupported region")


    if metadata:
        url = ''.join(['http://', host, ':', port,
                       "/AcuoREST/dicomrs/retrieve/studies/",
                        study])

        if series is not None:
            url += "/series/" + series
            if instance is not None:
                url += "/instances/" + instance

        url += "/metadata"

        r = requests.get(url, auth=(user, pw)) #, headers=headers

        if r.status_code != 200:
            print("Skipped series:", series)
            #raise ValueError("Invalid request (response code %d) for URL: %s" % (r.status_code, url))
            return None, url

        with open(filename, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

    else:
        url = ''.join(['http://', host, ':', port,
                       "/AcuoREST/wadoget?requestType=WADO&contentType=application/dicom&studyUID=",
                        study])

        if series is not None:
            url += "&seriesUID=" + series
            if instance is not None:
                url += "&objectUID=" + instance

        r = requests.get(url, auth=(user, pw)) #, headers=headers

        if r.status_code != 200:
            #print("Skipped instance:", instance)
            #raise ValueError("Invalid request (response code %d) for URL: %s" % (r.status_code, url))
            return None, url

        with open("temp.dcm", 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

        anonymize.anonymize(filename="temp.dcm", output_filename=filename)

        os.remove("temp.dcm")
            
        
    return r, url