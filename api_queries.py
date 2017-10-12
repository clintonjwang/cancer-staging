from dicom.examples import anonymize
import requests
import importlib
import os
import time

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


def retrieve_vna(user, pw, filename, study=None, series=None, instance=None, region='test', metadata=False, protocolExclude=None):
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

        url += "/metadata"+"?contentType=application/xml"

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


if __name__ == "__main__":
    user = ''
    pw = ''
    region = 'prod'

    accNum = 'TESTSIRIANO'
    r, url = search_vna(user, pw, region=region, accNum=accNum)
    study = r.json()[0]['0020000D']['Value'][0]

    r, url = search_vna(user, pw, region=region, study=study, modality=None)
    study_info = r.json()
    series = [ser['0020000E']['Value'][0] for ser in study_info]

    instances = {}

    for ser in series:
        r, url = search_vna(user, pw, region=region, study=study, series=ser)
        series_info = r.json()
        instances[ser] = [inst['00080018']['Value'][0] for inst in series_info]
        if len(instances[ser]) < 5:
            del instances[ser]

    base_dir = "raw_imgs\\" + accNum

    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
        
    total = 0
    skip_ser = 0
    skip_inst = 0
    for ser in instances:
        t = time.time()
        print("\n==============")
        print("Loading series", ser)
        img_dir = base_dir + "\\" + ser
        
        if not os.path.exists(img_dir):
            os.makedirs(img_dir)
            
        r, url = retrieve_vna(user, pw, region=region, filename = img_dir+"\\metadata.xml",
                          study=study, series=ser, metadata=True)
        if r is None:
            skip_ser += 1
            continue
            
        for count, inst in enumerate(instances[ser]):
            r, url = retrieve_vna(user, pw, region=region, filename = img_dir+"\\"+str(count)+".dcm",
                          study=study, series=ser, instance=inst)
            
            if r is not None:
                skip_inst += 1
            print(".", end="")
        
        total += count
        print("\nTime elapsed: %.2fs" % time.time()-t)
        
    print("Series loaded: ", len(series)-skip_ser, "/", len(series), sep="")
    print("Total images loaded:", total)
    print("Images skipped:", skip_inst)