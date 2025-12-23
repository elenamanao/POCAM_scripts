#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include "payload_header_drop.h"
#include "hdf5.h"
#define DS2482_800_MAX_CH 8
#define OW_UNUSED 0x03


#define MAX11614_SIPM2_CH 0
#define MAX11614_SIPM1_CH 1
#define MAX11614_MMB33_CH 2
#define MAX11614_MMB18_CH 3
#define MAX11614_KAPU2_CH 4
#define MAX11614_LMG2_CH  5
#define MAX11614_LMG1_CH  6
#define MAX11614_KAPU1_CH 7

static uint8_t *_buff;
static uint32_t n_left;

void printHeaderName(uint8_t rtype){
  char *name;
  switch (rtype){
    case PCM_START_OF_SUBRUN:
      name = "PCM_START_OF_SUBRUN";
      break;
    case PCM_CAL_TRIG_TIMES:
      name = "PCM_CAL_TRIG_TIMES";
      break;
    case PCM_BRD_TEMP_RECORD:
      name = "PCM_BRD_TEMP_RECORD";
      break;
    case PCM_PD_ADC_DATA:
      name = "PCM_PD_ADC_DATA";
      break;
    case PCM_TDC_DATA:
      name = "PCM_TDC_DATA";
      break;
    case PCM_END_OF_FILE:
      name = "PCM_END_OF_FILE";
      break;
    case PCM_DEVICEID_DATA:
      name ="PCM_DEVICEID_DATA";
      break;
    case PCM_FPGA_STATE_DATA:
      name ="PCM_FPGA_STATE_DATA";
      break;
    case PCM_VOLTAGE_DATA:
      name ="PCM_VOLTAGE_DATA";
      break;
    case PCM_MMB_SENSOR_DATA:
      name = "PCM_MMB_SENSOR_DATA";
      break;
    default:
      name = "";
      break;
  } 
  printf("Record type: 0x%02x (%s)\n",rtype, name  );
  return;
}

void printFPGAStateStream(uint8_t* buffer, uint32_t nbytes){
  printf("\tnbytes = %u; record size = %lu\n", nbytes ,  sizeof(pcm_fpga_state_t));
  if (nbytes != sizeof(pcm_fpga_state_t)){
    printf("\tWARNING! Size mismatch bytes\n");
  }
  pcm_fpga_state_t record;
  memcpy(&record, buffer, sizeof(pcm_fpga_state_t));
  char *name ;
  for (uint8_t i=0; i < 3; i++){
    switch (i){
      case 0:
        name = "IB";
        break;
      case 1:
        name = "DBM";
        break;
      case 2:
        name = "DBS";
        break;
      default:
        break;
    }
    printf("\t\t=========================\n");
    printf("\tTarget %u = %s\n", i, name);
    printf("\t\t------- registers -------\n");
    for (uint8_t j=0; j <16; j++){
      if (i==0){
        printf("\t\t reg=0x%x, val=0x%08x\n", j, record.ib_registers[j]);
      } else if (i==1){
        printf("\t\t reg=0x%x, val=0x%08x\n", j, record.dbm_registers[j]);
      } else if (i==2){
        printf("\t\t reg=0x%x, val=0x%08x\n", j, record.dbs_registers[j]);
      } 
    }
    printf("\t\t-------   pwms   --------\n");
    if (i == 0){
      for (uint8_t j = 0 ; j < 6 ; j++){
        printf("\t\t pwm=0x%x, val=%u\n", j, record.ib_pwms[j]);
      }
    } else if (i == 1){ 
      for (uint8_t j = 0 ; j < 4 ; j++){
        printf("\t\t pwm=0x%x, val=%u\n", j, record.dbm_pwms[j]);
      }
    } else if (i ==2 ) {
      for (uint8_t j = 0 ; j < 4 ; j++){
        printf("\t\t pwm=0x%x, val=%u\n", j, record.dbs_pwms[j]);
      }
    }
  }
  return;
}

void printDeviceIdStream(uint8_t* buffer, uint32_t nbytes){

  printf("\tnbytes = %u; record size = %lu\n", nbytes ,  sizeof(pcm_pocam_deviceid_t));
  if (nbytes != sizeof(pcm_pocam_deviceid_t)){
    printf("\tWARNING! Size mismatch bytes\n");
  }
  pcm_pocam_deviceid_t record;
  memcpy(&record, buffer, sizeof(pcm_pocam_deviceid_t));
  printf("\tMMB ID: 0x%s\n", record.mmb_id);
  printf("\tICM ID: 0x%s\n", record.icm_id);
  printf("\tICM FW: 0x%x\n", record.icm_fw_version);
  printf("\tIB  FPGA ID: 0x%016lx\n", record.pcm_fpga_traceids[0]);
  printf("\tDBM FPGA ID: 0x%016lx\n", record.pcm_fpga_traceids[1]);
  printf("\tDBS FPGA ID: 0x%016lx\n", record.pcm_fpga_traceids[2]);
  printf("\tIB  FPGA fw: 0x%x\n", record.pcm_fpga_fw_versions[0]);
  printf("\tDBM FPGA fw: 0x%x\n", record.pcm_fpga_fw_versions[1]);
  printf("\tDBS FPGA fw: 0x%x\n", record.pcm_fpga_fw_versions[2]);

  return;
}
void printTemperatureStream(uint8_t *buffer, uint32_t nbytes){
  uint16_t n_records = nbytes / sizeof(pcm_brd_temp_record_t);
  pcm_brd_temp_record_t record;
  for (uint16_t i = 0; i < n_records; i++){
    memcpy( &record, 
            buffer+i*sizeof(pcm_brd_temp_record_t), 
            sizeof(pcm_brd_temp_record_t));
    printf("\tch=%u, st=%u, raw=%04x, temp=%f\n", 
                record.channel, 
                record.status,
                record.raw_value,
                record.celcius
                );
  }
  uint32_t left = nbytes - n_records*sizeof(pcm_brd_temp_record_t);
  if ( left){
    printf("WARNING! %u extra bytes left (= %u - %ux%lu)\n", 
              left, nbytes, n_records, sizeof(pcm_brd_temp_record_t));
  }
  return;
}


void printSubrunStream(uint8_t *buffer, uint32_t nbytes){
  uint16_t n_records = nbytes / sizeof(pcm_subrun_start_t);
  if ( n_records != 1) {
    printf("Subrun stream should have only 1 record, has %u ... printing all \n", n_records);
  }
  pcm_subrun_start_t record;
  for (uint16_t i=0; i < n_records; i++){
    memcpy( &record, 
            buffer+i*sizeof(pcm_subrun_start_t), 
            sizeof(pcm_subrun_start_t));
    printf("\tsubrun=%u\n", record.subrun);
    printf("\ttrigs/subrun=%u, period=%u\n", record.trigs_per_subrun , record.trig_period);
    printf("\tts = %lu\n", record.start_time);
    printf("\tepoch ns = %lu\n", record.start_epoch_ns);
  }
  uint32_t left = nbytes - n_records*sizeof(pcm_subrun_start_t);
  if ( left){
    printf("WARNING! %u extra bytes left (= %u - %ux%lu)\n", 
              left, nbytes, n_records, sizeof(pcm_subrun_start_t));
  }
  return;
}

void printCalTrigStream(uint8_t *buffer, uint32_t nbytes){

  uint16_t n_records = nbytes / sizeof(uint64_t);
  uint64_t record;
  printf("N rec = %u\n", n_records);
  for (uint16_t i=0; i < n_records; i++){
    memcpy( &record, 
            buffer+i*sizeof(uint64_t), 
            sizeof(uint64_t));
    printf("\ticm ts=%lu\n",record);
  }

}

void printPDADCStream(uint8_t *buffer, uint32_t nbytes){
  uint16_t n_records = (nbytes - sizeof(pcm_pd_adc_header_t)) / 
                                      sizeof(pcm_pd_adc_record_t);
  pcm_pd_adc_header_t header;
  memcpy(&header, buffer, sizeof(pcm_pd_adc_header_t));
  printf("source=0x%02x, n=%u, avg=%u, w=%u \n", 
                                header.source,
                                header.n_records, 
                                header.averaging, 
                                header.window);
  pcm_pd_adc_record_t record;
  for (uint16_t i=0; i < n_records; i++){
    memcpy(&record, 
           buffer + sizeof(pcm_pd_adc_header_t) + i*sizeof(pcm_pd_adc_record_t),
           sizeof(pcm_pd_adc_record_t));
    printf("\t base=%u, max=%u, area=%u \n", record.base, record.max_ampl, record.area);
  
  }
  return;

}

void printTDCStream(uint8_t *buffer, uint32_t nbytes){
  uint16_t n_records = (nbytes - sizeof(pcm_tdc_header_t)) / 
                                 sizeof(uint64_t);
  pcm_tdc_header_t header;
  memcpy(&header, buffer, sizeof(pcm_tdc_header_t));
  printf("source=0x%02x, n=%u, m=0x%x, hv=%u, thr=%u \n", 
                                header.source,
                                header.n_records, 
                                header.mode, 
                                header.hv,
                                header.threshold
                                );                             
  pcm_tdc_record_t record;
  double prev_time = 0.0;
  double cur_time=0.0;
  double epoch_delay = 0.0;
  for (uint16_t i=0; i < n_records; i++){
    memcpy(&record, 
           buffer + sizeof(pcm_tdc_header_t) + i*sizeof(pcm_tdc_record_t),
           sizeof(pcm_tdc_record_t));
           
    cur_time=record.coarse/0.175 + record.fine/(0.350*16.);
    printf("\tch=%u  ro=0x%02x e=%u c=%u f=%u" , 
        record.channel, 
        record.ro_err_bits, 
        record.edge, 
        record.coarse, 
        record.fine
       );

    if (record.ro_err_bits == 0x10){
      epoch_delay+=((1<<29) / 0.175);
      printf("\n");
    } else {
      printf(" dt = %0.3f\n",  epoch_delay+cur_time - prev_time);
      epoch_delay=0.0;
      prev_time = cur_time;
    }
  } 
  return;
}

void printMMBSensorSteam(uint8_t *buffer, uint32_t nbytes){
  if (nbytes < sizeof(pcm_mmb_sensor_record_t)){
    printf("\t WARNING! Not enough bytes in record! Skipping...\n");
  }
  if (nbytes >  sizeof(pcm_mmb_sensor_record_t)){
    printf("Too many bytes in record! Caution adviced\n");
  }
  pcm_mmb_sensor_record_t record;
  memcpy(&record, buffer, sizeof(pcm_mmb_sensor_record_t));
  printf("\t Accel. = [ %0.3f, %0.3f, %0.3f ] m/s2, T = %0.2fC\n",
                                    record.accelerometerXYZ[0], 
                                    record.accelerometerXYZ[1], 
                                    record.accelerometerXYZ[2],
                                    record.accelerometerT);
  printf("\t MagF. = [ %0.3f, %0.3f, %0.3f ] muT, T = %0.2fC\n",
                                    record.magnetometerXYZ[0]*1e6, 
                                    record.magnetometerXYZ[1]*1e6, 
                                    record.magnetometerXYZ[2]*1e6,
                                    record.magnetometerT);
  printf("\t P = %0.3f hPa, T = %0.2fC\n",
                                    record.pressureP,
                                    record.pressureT);
  return;
}

void printVoltageStream(uint8_t *buffer, uint32_t nbytes){
  pcm_voltage_record_t record;
  uint32_t n_rec = nbytes/sizeof(pcm_voltage_record_t);
  printf("\tn recods = %u\n", n_rec);
  if ( (nbytes - n_rec * sizeof(pcm_voltage_record_t))>0){
    printf("WARNING! Extra bytes left!");
  }
  char *name;
  for (uint8_t i =0 ; i< n_rec; i++){
    memcpy(&record, 
           buffer + i*sizeof(pcm_voltage_record_t), 
           sizeof(pcm_voltage_record_t));
    switch( record.channel){
      case MAX11614_MMB18_CH:
        name = "mmb18";
        break;
      case MAX11614_MMB33_CH:
        name = "mmb33";
        break;
      case MAX11614_KAPU1_CH:
        name = "kapu1";
        break;
      case MAX11614_KAPU2_CH:
        name = "kapu2";
        break;
      case MAX11614_LMG1_CH:
        name = "lmg1";
        break;
      case MAX11614_LMG2_CH:
        name = "lmg2";
        break;
      case MAX11614_SIPM1_CH:
        name = "sipm1";
        break;
      case MAX11614_SIPM2_CH:
        name = "sipm2";
        break;
      default:
        name="";
        break; 
    }
    printf("\t\t name=%s\tV = %0.3f V\t(ch=%u, raw=%u)\n", 
            name,
            record.voltage,
            record.channel, 
            record.raw);
  }
  return;
}

void printRecord(common_header_t header, FILE *fptr){
  printHeaderName(header.record_type);
  uint64_t icm_ts = header.icm_ts_low;
  icm_ts |= (((uint64_t)(header.icm_ts_high)) << 32); 
  printf("ICM ts:  %lu \n", icm_ts);
  
  uint32_t nbytes = header.nbytes_low;
  nbytes |= (((uint32_t)(header.nbytes_high)) << 16);
  printf("Nbytes: %u \n", nbytes);
  
  uint32_t nb_record = nbytes -sizeof(common_header_t);
  if (nb_record==0) { return; }
  if ( header.record_type == PCM_END_OF_FILE){
    printf("End of file record: finishing");
    n_left = 0;
    return;
  }
  if (nb_record > n_left){
    printf("Erorr! Requested %u , available %u bytes \n",nb_record , n_left);
    return;
  }
  fread(_buff, nb_record,1,fptr);
  n_left = n_left - nb_record;
  
  switch (header.record_type){
    case PCM_START_OF_SUBRUN :
      printSubrunStream(_buff, nb_record);
      break;
    case PCM_BRD_TEMP_RECORD :
      printTemperatureStream(_buff, nb_record);
      break;
    case PCM_CAL_TRIG_TIMES :
      printCalTrigStream(_buff, nb_record);
      break;
    case PCM_PD_ADC_DATA:
      printPDADCStream(_buff, nb_record);
      break;
    case PCM_TDC_DATA:
      printTDCStream(_buff, nb_record);
      break;
    case PCM_DEVICEID_DATA:
      printDeviceIdStream(_buff, nb_record);
      break;
    case PCM_FPGA_STATE_DATA:
      printFPGAStateStream(_buff, nb_record);
      break;
    case PCM_MMB_SENSOR_DATA:
      printMMBSensorSteam(_buff, nb_record);
      break;
    case PCM_VOLTAGE_DATA:
      printVoltageStream(_buff, nb_record);
      break;
    default:
      break;
  }
  
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s file1.pdd [file2.pdd ...]\n", argv[0]);
        return 1;
    }

    _buff = calloc(1, sizeof(uint64_t) * 2048 + 1024);
    if (!_buff) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    // --- Create a single HDF5 output file ---
    hid_t file_id = H5Fcreate("combined_subruns.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error creating combined_subruns.h5\n");
        free(_buff);
        return 1;
    }

    // --- Process each input file ---
    for (int fidx = 1; fidx < argc; fidx++) {
        char *fname = argv[fidx];
        FILE *ptr = fopen(fname, "rb");
        if (!ptr) {
            fprintf(stderr, "Cannot open %s\n", fname);
            continue;
        }

        printf("Processing %s...\n", fname);

        // Create subrun group for this file (iterating by index)
        char gname[64];
        snprintf(gname, sizeof(gname), "subrun_%03d", fidx - 1);
        hid_t grp = H5Gcreate2(file_id, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    fseek(ptr, 0, SEEK_END);
    n_left = ftell(ptr);
    fseek(ptr, 0, SEEK_SET);

    // Counters for unique dataset/group names
    int subrun_count = 0, caltrig_count = 0, temp_count = 0, adc_count = 0;
    int tdc_count = 0, devid_count = 0, fpga_count = 0, mmb_count = 0, volt_count = 0;

    common_header_t header;
    while (n_left >= sizeof(common_header_t)) {
        fread(&header, sizeof(common_header_t), 1, ptr);
        n_left -= sizeof(common_header_t);

        uint32_t nbytes = header.nbytes_low | ((uint32_t)header.nbytes_high << 16);
        if (nbytes <= sizeof(common_header_t)) continue;
        fread(_buff, nbytes - sizeof(common_header_t), 1, ptr);
        n_left -= (nbytes - sizeof(common_header_t));

        switch (header.record_type) {

            // ---------- Subrun ----------
            case PCM_START_OF_SUBRUN: {
                char dname[64];
                snprintf(dname, sizeof(dname), "subrun_%03d", subrun_count++);
                hsize_t dims[1] = {1};
                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_subrun_start_t));
                H5Tinsert(dtype, "subrun", HOFFSET(pcm_subrun_start_t, subrun), H5T_NATIVE_UINT16);
                H5Tinsert(dtype, "trigs_per_subrun", HOFFSET(pcm_subrun_start_t, trigs_per_subrun), H5T_NATIVE_UINT16);
                H5Tinsert(dtype, "trig_period", HOFFSET(pcm_subrun_start_t, trig_period), H5T_NATIVE_UINT16);
                H5Tinsert(dtype, "start_time", HOFFSET(pcm_subrun_start_t, start_time), H5T_NATIVE_UINT64);
                H5Tinsert(dtype, "start_epoch_ns", HOFFSET(pcm_subrun_start_t, start_epoch_ns), H5T_NATIVE_UINT64);
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dset = H5Dcreate2(grp, dname, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, _buff);
                H5Dclose(dset); H5Sclose(space); H5Tclose(dtype);
                break;
            }

            // ---------- Cal Trig Times ----------
            case PCM_CAL_TRIG_TIMES: {
                uint16_t nrec = (nbytes - sizeof(common_header_t)) / sizeof(uint64_t);
                hsize_t dims[1] = { nrec };
                char dname[64];
                snprintf(dname, sizeof(dname), "cal_trig_%03d", caltrig_count++);
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dset = H5Dcreate2(grp, dname, H5T_NATIVE_UINT64, space,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, _buff);
                H5Dclose(dset); H5Sclose(space);
                break;
            }

            // ---------- Temperature ----------
            case PCM_BRD_TEMP_RECORD: {
                uint16_t nrec = (nbytes - sizeof(common_header_t)) / sizeof(pcm_brd_temp_record_t);
                hsize_t dims[1] = { nrec };
                char dname[64];
                snprintf(dname, sizeof(dname), "temperature_%03d", temp_count++);
                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_brd_temp_record_t));
                H5Tinsert(dtype, "channel", HOFFSET(pcm_brd_temp_record_t, channel), H5T_NATIVE_UINT8);
                H5Tinsert(dtype, "status",  HOFFSET(pcm_brd_temp_record_t, status),  H5T_NATIVE_UINT8);
                H5Tinsert(dtype, "raw_value", HOFFSET(pcm_brd_temp_record_t, raw_value), H5T_NATIVE_UINT16);
                H5Tinsert(dtype, "celcius", HOFFSET(pcm_brd_temp_record_t, celcius), H5T_NATIVE_FLOAT);
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dset = H5Dcreate2(grp, dname, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, _buff);
                H5Dclose(dset); H5Sclose(space); H5Tclose(dtype);
                break;
            }

            // ---------- PD ADC ----------
            case PCM_PD_ADC_DATA: {
                char gname[64];
                snprintf(gname, sizeof(gname), "adc_%03d", adc_count++);
                hid_t grps = H5Gcreate2(grp, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                pcm_pd_adc_header_t header_adc;
                memcpy(&header_adc, _buff, sizeof(pcm_pd_adc_header_t));
                hid_t htype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_pd_adc_header_t));
                H5Tinsert(htype, "source", HOFFSET(pcm_pd_adc_header_t, source), H5T_NATIVE_UINT8);
                H5Tinsert(htype, "n_records", HOFFSET(pcm_pd_adc_header_t, n_records), H5T_NATIVE_UINT16);
                H5Tinsert(htype, "averaging", HOFFSET(pcm_pd_adc_header_t, averaging), H5T_NATIVE_UINT8);
                H5Tinsert(htype, "window", HOFFSET(pcm_pd_adc_header_t, window), H5T_NATIVE_UINT16);
                hsize_t hdims[1] = {1};
                hid_t hspace = H5Screate_simple(1, hdims, NULL);
                hid_t hset = H5Dcreate2(grps, "header", htype, hspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(hset, htype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header_adc);
                H5Dclose(hset); H5Sclose(hspace); H5Tclose(htype);

                uint16_t nrec = (nbytes - sizeof(common_header_t) - sizeof(pcm_pd_adc_header_t)) / sizeof(pcm_pd_adc_record_t);
                hsize_t dims[1] = { nrec };
                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_pd_adc_record_t));
                H5Tinsert(dtype, "area", HOFFSET(pcm_pd_adc_record_t, area), H5T_NATIVE_UINT32);
                H5Tinsert(dtype, "max_ampl", HOFFSET(pcm_pd_adc_record_t, max_ampl), H5T_NATIVE_UINT16);
                H5Tinsert(dtype, "base", HOFFSET(pcm_pd_adc_record_t, base), H5T_NATIVE_UINT16);
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dset = H5Dcreate2(grps, "records", dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, _buff + sizeof(pcm_pd_adc_header_t));
                H5Dclose(dset); H5Sclose(space); H5Tclose(dtype);
                H5Gclose(grps);
                break;
            }

            // ---------- TDC ----------
            case PCM_TDC_DATA: {
                char gname[64];
                snprintf(gname, sizeof(gname), "tdc_%03d", tdc_count++);
                hid_t grps = H5Gcreate2(grp, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                pcm_tdc_header_t header_tdc;
                memcpy(&header_tdc, _buff, sizeof(pcm_tdc_header_t));
                hid_t htype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_tdc_header_t));
                H5Tinsert(htype, "source", HOFFSET(pcm_tdc_header_t, source), H5T_NATIVE_UINT8);
                H5Tinsert(htype, "n_records", HOFFSET(pcm_tdc_header_t, n_records), H5T_NATIVE_UINT16);
                H5Tinsert(htype, "mode", HOFFSET(pcm_tdc_header_t, mode), H5T_NATIVE_UINT8);
                H5Tinsert(htype, "hv", HOFFSET(pcm_tdc_header_t, hv), H5T_NATIVE_UINT8);
                H5Tinsert(htype, "threshold", HOFFSET(pcm_tdc_header_t, threshold), H5T_NATIVE_UINT16);
                hsize_t hdims[1] = {1};
                hid_t hspace = H5Screate_simple(1, hdims, NULL);
                hid_t hset = H5Dcreate2(grps, "header", htype, hspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(hset, htype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header_tdc);
                H5Dclose(hset); H5Sclose(hspace); H5Tclose(htype);

                uint16_t nrec = (nbytes - sizeof(common_header_t) - sizeof(pcm_tdc_header_t)) / sizeof(pcm_tdc_record_t);
                hsize_t dims[1] = { nrec };
                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_tdc_record_t));
                H5Tinsert(dtype, "channel", HOFFSET(pcm_tdc_record_t, channel), H5T_NATIVE_UINT8);
                H5Tinsert(dtype, "ro_err_bits", HOFFSET(pcm_tdc_record_t, ro_err_bits), H5T_NATIVE_UINT8);
                H5Tinsert(dtype, "edge", HOFFSET(pcm_tdc_record_t, edge), H5T_NATIVE_UINT8);
                H5Tinsert(dtype, "coarse", HOFFSET(pcm_tdc_record_t, coarse), H5T_NATIVE_UINT32);
                H5Tinsert(dtype, "fine", HOFFSET(pcm_tdc_record_t, fine), H5T_NATIVE_UINT8);
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dset = H5Dcreate2(grps, "records", dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, _buff + sizeof(pcm_tdc_header_t));
                H5Dclose(dset); H5Sclose(space); H5Tclose(dtype);
                H5Gclose(grps);
                break;
            }

            // ---------- Device ID ----------
            case PCM_DEVICEID_DATA: {
                char dname[64];
                snprintf(dname, sizeof(dname), "deviceid_%03d", devid_count++);
                hsize_t dims[1] = {1};
                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_pocam_deviceid_t));
                H5Tinsert(dtype, "icm_id", HOFFSET(pcm_pocam_deviceid_t, icm_id), H5T_C_S1);
                H5Tinsert(dtype, "mmb_id", HOFFSET(pcm_pocam_deviceid_t, mmb_id), H5T_C_S1);
                H5Tinsert(dtype, "icm_fw_version", HOFFSET(pcm_pocam_deviceid_t, icm_fw_version), H5T_NATIVE_UINT16);
                H5Tinsert(dtype, "pcm_fpga_traceids", HOFFSET(pcm_pocam_deviceid_t, pcm_fpga_traceids), H5Tarray_create2(H5T_NATIVE_UINT64, 1, (hsize_t[]){3}));
                H5Tinsert(dtype, "pcm_fpga_fw_versions", HOFFSET(pcm_pocam_deviceid_t, pcm_fpga_fw_versions), H5Tarray_create2(H5T_NATIVE_UINT16, 1, (hsize_t[]){3}));
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dset = H5Dcreate2(grp, dname, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, _buff);
                H5Dclose(dset); H5Sclose(space); H5Tclose(dtype);
                break;
            }

            // ---------- FPGA State ----------
            case PCM_FPGA_STATE_DATA: {
                char dname[64];
                snprintf(dname, sizeof(dname), "fpga_state_%03d", fpga_count++);
                hsize_t dims[1] = {1};
                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_fpga_state_t));
                H5Tinsert(dtype, "ib_registers", HOFFSET(pcm_fpga_state_t, ib_registers), H5Tarray_create2(H5T_NATIVE_UINT32, 1, (hsize_t[]){16}));
                H5Tinsert(dtype, "dbm_registers", HOFFSET(pcm_fpga_state_t, dbm_registers), H5Tarray_create2(H5T_NATIVE_UINT32, 1, (hsize_t[]){16}));
                H5Tinsert(dtype, "dbs_registers", HOFFSET(pcm_fpga_state_t, dbs_registers), H5Tarray_create2(H5T_NATIVE_UINT32, 1, (hsize_t[]){16}));
                H5Tinsert(dtype, "ib_pwms", HOFFSET(pcm_fpga_state_t, ib_pwms), H5Tarray_create2(H5T_NATIVE_UINT32, 1, (hsize_t[]){6}));
                H5Tinsert(dtype, "dbm_pwms", HOFFSET(pcm_fpga_state_t, dbm_pwms), H5Tarray_create2(H5T_NATIVE_UINT32, 1, (hsize_t[]){4}));
                H5Tinsert(dtype, "dbs_pwms", HOFFSET(pcm_fpga_state_t, dbs_pwms), H5Tarray_create2(H5T_NATIVE_UINT32, 1, (hsize_t[]){4}));
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dset = H5Dcreate2(grp, dname, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, _buff);
                H5Dclose(dset); H5Sclose(space); H5Tclose(dtype);
                break;
            }

            // ---------- MMB Sensor ----------
            case PCM_MMB_SENSOR_DATA: {
                char dname[64];
                snprintf(dname, sizeof(dname), "mmb_sensor_%03d", mmb_count++);
                hsize_t dims[1] = {1};
                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_mmb_sensor_record_t));
                H5Tinsert(dtype, "accelerometerXYZ", HOFFSET(pcm_mmb_sensor_record_t, accelerometerXYZ), H5Tarray_create2(H5T_NATIVE_FLOAT, 1, (hsize_t[]){3}));
                H5Tinsert(dtype, "accelerometerT", HOFFSET(pcm_mmb_sensor_record_t, accelerometerT), H5T_NATIVE_FLOAT);
                H5Tinsert(dtype, "magnetometerXYZ", HOFFSET(pcm_mmb_sensor_record_t, magnetometerXYZ), H5Tarray_create2(H5T_NATIVE_FLOAT, 1, (hsize_t[]){3}));
                H5Tinsert(dtype, "magnetometerT", HOFFSET(pcm_mmb_sensor_record_t, magnetometerT), H5T_NATIVE_FLOAT);
                H5Tinsert(dtype, "pressureP", HOFFSET(pcm_mmb_sensor_record_t, pressureP), H5T_NATIVE_FLOAT);
                H5Tinsert(dtype, "pressureT", HOFFSET(pcm_mmb_sensor_record_t, pressureT), H5T_NATIVE_FLOAT);
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dset = H5Dcreate2(grp, dname, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, _buff);
                H5Dclose(dset); H5Sclose(space); H5Tclose(dtype);
                break;
            }

            // ---------- Voltage ----------
            case PCM_VOLTAGE_DATA: {
                uint16_t nrec = (nbytes - sizeof(common_header_t)) / sizeof(pcm_voltage_record_t);
                hsize_t dims[1] = { nrec };
                char dname[64];
                snprintf(dname, sizeof(dname), "voltage_%03d", volt_count++);
                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(pcm_voltage_record_t));
                H5Tinsert(dtype, "channel", HOFFSET(pcm_voltage_record_t, channel), H5T_NATIVE_UINT8);
                H5Tinsert(dtype, "raw", HOFFSET(pcm_voltage_record_t, raw), H5T_NATIVE_UINT16);
                H5Tinsert(dtype, "voltage", HOFFSET(pcm_voltage_record_t, voltage), H5T_NATIVE_FLOAT);
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dset = H5Dcreate2(grp, dname, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, _buff);
                H5Dclose(dset); H5Sclose(space); H5Tclose(dtype);
                break;
            }

            default:
                break;
        }
    }

    H5Gclose(grp);
    fclose(ptr);
    printf("Finished %s -> %s\n", fname, gname);
}

H5Fclose(file_id);
free(_buff);
printf("All subruns saved to combined_subruns.h5\n");
return 0;
}




